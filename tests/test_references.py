import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

import _gseapy_pre_rank_wrap as wrapper

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('reference_downloader', ROOT / 'config/download_gseapy_gmt_files.py')
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def mapping(self, rows, target='mouse', name='mapping.csv'):
        path = self.root / name
        pd.DataFrame(rows, columns=['gene_name', target + '_gene_name', target + '_gene_names_all']).to_csv(path, index=False)
        return path

    def test_strict_mapping_qc_and_empty_omission(self):
        mapping = self.mapping([('A','a','a'), ('A','a','a'), ('B','b','b, c'),
                                ('C','d','d'), ('D','d','d'), ('E','','')])
        source = self.root / 'source.gmt'; source.write_text('T\tdesc\tA\tA\tB\tC\tE\tunknown\nEmpty\tx\tB\n')
        output = self.root / 'converted.gmt'
        wrapper.convert_gmt_homologs(str(source), str(output), homolog_mapping_h2m_file=str(mapping))
        self.assertEqual(output.read_text(), 'T\tdesc\ta\n')
        qc = json.loads(Path(str(output)+'.qc.json').read_text())
        self.assertEqual(qc['totals'], dict(source_count=7, mapped_count=2, excluded_count=3,
                                          unmapped_count=2, deduplicated_count=1, output_count=1,
                                          source_term_count=2, output_term_count=1, empty_term_count=1))
        self.assertEqual(qc['mapping_policy']['one_to_many_source_count'], 1)
        self.assertEqual(qc['mapping_policy']['many_to_one_source_count'], 2)
        self.assertEqual(qc['mapping_policy']['duplicate_relationship_record_count'], 1)
        self.assertEqual(qc['source_sha256'], wrapper._sha256_file(source))
        self.assertEqual(qc['output_sha256'], wrapper._sha256_file(output))

    def test_mouse_to_human_and_complete_columns_required(self):
        mapping = self.mapping([('a','A','A')], target='human')
        source = self.root / 'in.gmt'; source.write_text('T\t\ta\n')
        output = self.root / 'out.gmt'
        wrapper.convert_gmt_homologs(str(source), str(output), covert2organism='Human', homolog_mapping_m2h_file=str(mapping))
        self.assertEqual(output.read_text(), 'T\t\tA\n')
        mapping.write_text('gene_name,human_gene_name\na,A\n')
        with self.assertRaisesRegex(ValueError, 'complete relationship'):
            wrapper.convert_gmt_homologs(str(source), str(self.root/'new.gmt'), covert2organism='Human', homolog_mapping_m2h_file=str(mapping))
        self.assertFalse((self.root/'new.gmt').exists())

    def test_conflicting_source_records_reject_entire_map_before_outputs(self):
        mapping = self.mapping([('A','a','a'), ('Z','z','z'), ('Z','x','x')])
        source = self.root/'in.gmt'; source.write_text('T\t\tA\n')
        output = self.root/'new'/'out.gmt'
        with self.assertRaisesRegex(ValueError, 'Contradictory.*Z'):
            wrapper.convert_gmt_homologs(str(source), str(output), homolog_mapping_h2m_file=str(mapping))
        self.assertFalse(output.parent.exists())

    def declaration(self, species=('mouse','human')):
        mapping = self.mapping([('A','a','a')])
        resources = [dict(id='resource'+str(index), provider='enrichr', collection='library'+str(index),
                          identifier_species=organism, identifier_species_evidence_url='https://example.org/evidence',
                          source_url='https://example.org/source', endpoint_organism='Human', release=None)
                     for index,organism in enumerate(species)]
        data = dict(schema_version=1, target_species='mouse', resources=resources, mappings={'human_to_mouse':mapping.name})
        path=self.root/'declaration.json'; path.write_text(json.dumps(data))
        return path, data

    def fake_download(self, requested, outdir, organism):
        path=Path(outdir)/(requested[0]+'.gmt'); path.parent.mkdir(parents=True)
        # Endpoint is deliberately identical while identifier declarations differ.
        path.write_text('T\tdesc\t'+('a' if requested[0]=='library0' else 'A')+'\n')
        return [str(path)]

    def test_native_bytes_conversion_and_deterministic_manifest(self):
        declaration,_=self.declaration()
        with patch.object(downloader,'make_gmts_from_enrichr_libraries',side_effect=self.fake_download):
            first=downloader.acquire_snapshot(declaration,self.root/'one')
            second=downloader.acquire_snapshot(declaration,self.root/'two')
        self.assertEqual((first/'manifest.json').read_bytes(),(second/'manifest.json').read_bytes())
        self.assertEqual((first/'source/resource0/library0.gmt').read_bytes(),(first/'routed/mouse/resource0.gmt').read_bytes())
        self.assertEqual((first/'routed/mouse/resource1.gmt').read_text(),'T\tdesc\ta\n')
        manifest=json.loads((first/'manifest.json').read_text())
        self.assertEqual([item['route'] for item in manifest['resources']],['native','human_to_mouse'])
        self.assertIsNone(manifest['resources'][1]['mapping_release'])
        self.assertEqual(manifest['resources'][1]['mapping_sha256'],wrapper._sha256_file(self.root/'mapping.csv'))

    def test_unknown_species_blocks_all_before_download(self):
        declaration,data=self.declaration()
        data['resources'][1]['identifier_species']='mixed'
        declaration.write_text(json.dumps(data))
        with patch.object(downloader,'make_gmts_from_enrichr_libraries') as download:
            with self.assertRaisesRegex(ValueError,'unknown/mixed'):
                downloader.acquire_snapshot(declaration,self.root/'output')
        download.assert_not_called()
        self.assertFalse((self.root/'output').exists())

    def test_interruption_and_existing_snapshot_preserved(self):
        declaration,_=self.declaration()
        old=self.root/'old'; old.mkdir(); (old/'sentinel').write_text('original')
        with patch.object(downloader,'make_gmts_from_enrichr_libraries') as download:
            with self.assertRaises(FileExistsError):
                downloader.acquire_snapshot(declaration,old)
        download.assert_not_called()
        for failure in (RuntimeError('download failed'), KeyboardInterrupt()):
            calls=[]
            def interrupted(*args,**kwargs):
                calls.append(1)
                if len(calls)>1:
                    raise failure
                return self.fake_download(*args,**kwargs)
            with patch.object(downloader,'make_gmts_from_enrichr_libraries',side_effect=interrupted):
                with self.assertRaises(type(failure)):
                    downloader.acquire_snapshot(declaration,self.root/'new')
            self.assertFalse((self.root/'new').exists())
            self.assertFalse(list(self.root.glob('.new.staging-*')))
        self.assertEqual((old/'sentinel').read_text(),'original')

    def test_conflicting_scalar_with_same_complete_relations_rejected(self):
        mapping=self.mapping([('A','a','a, b'),('A','b','b, a')])
        with self.assertRaisesRegex(ValueError,'Contradictory'):
            wrapper._load_homolog_mapping(str(mapping),'gene_name','mouse_gene_name')

    def test_publish_failure_cleans_staging(self):
        declaration,_=self.declaration()
        with patch.object(downloader,'make_gmts_from_enrichr_libraries',side_effect=self.fake_download), patch.object(downloader.os,'rename',side_effect=OSError('interrupted publication')):
            with self.assertRaises(OSError):
                downloader.acquire_snapshot(declaration,self.root/'new')
        self.assertFalse((self.root/'new').exists())
        self.assertFalse(list(self.root.glob('.new.staging-*')))

    def test_msigdb_requested_missing_not_silently_omitted(self):
        with patch.object(wrapper.gmsigdb.Msigdb,'get_gmt',return_value=None):
            with self.assertRaisesRegex(ValueError,'requested category'):
                wrapper.make_gmts_from_Msigdb('2025.1.Hs',['missing'],str(self.root/'download'))

    def test_msigdb_acquisition_human_target_and_release_gate(self):
        declaration,data=self.declaration(species=('human','mouse'))
        mapping=self.mapping([('a','A','A')],target='human',name='reverse.csv')
        data['target_species']='human'
        data['mappings']={'mouse_to_human':mapping.name}
        for resource,release in zip(data['resources'],('2025.1.Hs','2025.1.Mm')):
            resource.update(provider='msigdb',release=release)
        declaration.write_text(json.dumps(data))
        def backend(category,dbver,entrez):
            self.assertFalse(entrez)
            return {'T':['A' if dbver.endswith('Hs') else 'a']}
        with patch.object(wrapper.gmsigdb.Msigdb,'get_gmt',side_effect=backend) as fetch:
            output=downloader.acquire_snapshot(declaration,self.root/'human')
        self.assertEqual([call.kwargs['dbver'] for call in fetch.call_args_list],['2025.1.Hs','2025.1.Mm'])
        self.assertEqual([call.kwargs['category'] for call in fetch.call_args_list],['library0','library1'])
        self.assertEqual((output/'source/resource0/library0.gmt').read_bytes(),(output/'routed/human/resource0.gmt').read_bytes())
        self.assertTrue((output/'routed/human/resource1.gmt').read_text().endswith('\tA\n'))
        data['resources'][1]['release']='2025.1.Hs'
        declaration.write_text(json.dumps(data))
        with patch.object(wrapper.gmsigdb.Msigdb,'get_gmt') as fetch:
            with self.assertRaisesRegex(ValueError,'expected approved release'):
                downloader.acquire_snapshot(declaration,self.root/'invalid')
        fetch.assert_not_called()


if __name__=='__main__':
    unittest.main()
