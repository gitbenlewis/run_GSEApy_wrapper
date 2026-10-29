#!/usr/bin/env python3
"""Acquire a declared reference snapshot into a new, atomically published directory."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _gseapy_pre_rank_wrap import (
    make_gmts_from_enrichr_libraries, make_gmts_from_Msigdb, convert_gmt_homologs,
    _homolog_policy, _read_gmt_records, _sha256_file,
)


def acquire_snapshot(manifest_path, output):
    """Publish all declared resources or none; never replace an old snapshot."""
    manifest_path, output = Path(manifest_path).resolve(), Path(output).absolute()
    declaration = json.loads(manifest_path.read_text())
    if declaration.get('schema_version') != 1:
        raise ValueError('Reference declaration schema_version must be 1.')
    target = declaration.get('target_species')
    if target not in ('human', 'mouse'):
        raise ValueError('target_species must be human or mouse.')
    resources = declaration.get('resources')
    if not isinstance(resources, list) or not resources:
        raise ValueError('Declare every requested resource in a nonempty resources list.')
    if output.exists() or output.is_symlink():
        raise FileExistsError(f'Refusing to replace snapshot {output}; choose a new output directory.')
    identifiers, mappings = set(), {}
    for resource in resources:
        identifier = resource.get('id', '')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', identifier) or identifier in identifiers:
            raise ValueError(f'Resource id must be unique and filename-safe: {identifier!r}')
        identifiers.add(identifier)
        species = resource.get('identifier_species')
        if species not in ('human', 'mouse'):
            raise ValueError(f'{identifier}: unknown/mixed identifier species blocks routing; provide an evidence-backed declaration.')
        for key in ('identifier_species_evidence_url', 'source_url'):
            if not isinstance(resource.get(key), str) or not resource[key].startswith(('https://', 'http://')):
                raise ValueError(f'{identifier}: provide {key}; evidence is caller-declared, not automatically verified.')
        provider = resource.get('provider')
        if provider not in ('enrichr', 'msigdb'):
            raise ValueError(f'{identifier}: provider must be enrichr or msigdb.')
        collection = resource.get('collection', '')
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', collection):
            raise ValueError(f'{identifier}: collection must be an explicit filename-safe provider key.')
        if provider == 'enrichr' and resource.get('endpoint_organism') not in ('Human', 'Mouse'):
            raise ValueError(f'{identifier}: declare endpoint_organism independently of identifier species.')
        if provider == 'msigdb':
            expected = '2025.1.Hs' if species == 'human' else '2025.1.Mm'
            if resource.get('release') != expected:
                raise ValueError(f'{identifier}: expected approved release {expected} for declared identifiers.')
        if species != target:
            direction = species + '_to_' + target
            path = declaration.get('mappings', {}).get(direction)
            if not path:
                raise ValueError(f'{identifier}: provide mappings.{direction} with complete relationship columns.')
            path = Path(path)
            path = path if path.is_absolute() else manifest_path.parent / path
            _homolog_policy(path, 'gene_name', 'mouse_gene_name' if target == 'mouse' else 'human_gene_name')
            mappings[direction] = path
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='.' + output.name + '.staging-', dir=output.parent))
    try:
        records = []
        for resource in resources:
            identifier, species = resource['id'], resource['identifier_species']
            source_dir = stage / 'source' / identifier
            if resource['provider'] == 'enrichr':
                files = make_gmts_from_enrichr_libraries([resource['collection']], str(source_dir), organism=resource['endpoint_organism'])
            else:
                files = make_gmts_from_Msigdb(resource['release'], [resource['collection']], str(source_dir))
            expected = source_dir / (resource['collection'] + '.gmt')
            if len(files) != 1 or Path(files[0]) != expected or not expected.is_file():
                raise ValueError(f'{identifier}: requested download did not produce exactly its declared GMT.')
            source_records = _read_gmt_records(expected)
            if not source_records:
                raise ValueError(f'{identifier}: downloaded GMT has no terms.')
            routed = stage / 'routed' / target / (identifier + '.gmt')
            routed.parent.mkdir(parents=True, exist_ok=True)
            mapping_hash = None
            if species == target:
                shutil.copyfile(expected, routed)
                qc = dict(native_preserved=True,
                          terms=[dict(term=term, source_count=len(genes), output_count=len(genes),
                                      source_unique_count=len(set(genes)), empty=False)
                                 for term, _, genes in source_records],
                          totals=dict(source_term_count=len(source_records), output_term_count=len(source_records),
                                      source_count=sum(len(genes) for _, _, genes in source_records),
                                      output_count=sum(len(genes) for _, _, genes in source_records)))
            else:
                mapping = mappings[species + '_to_' + target]
                options = {'homolog_mapping_h2m_file' if target == 'mouse' else 'homolog_mapping_m2h_file': str(mapping)}
                convert_gmt_homologs(str(expected), str(routed), covert2organism=target, **options)
                qc = json.loads(Path(str(routed) + '.qc.json').read_text())
                mapping_hash = _sha256_file(mapping)
            record = {key: resource.get(key) for key in
                      ('id', 'provider', 'collection', 'release', 'endpoint_organism', 'source_url', 'source_study_species',
                       'identifier_species', 'identifier_species_evidence_url', 'license', 'license_url')}
            record.update(species_evidence_status='caller_declared_not_automatically_verified',
                          route='native' if species == target else species + '_to_' + target,
                          target_species=target, source_path=str(expected.relative_to(stage)),
                          output_path=str(routed.relative_to(stage)), source_sha256=_sha256_file(expected),
                          output_sha256=_sha256_file(routed), mapping_sha256=mapping_hash,
                          mapping_source_url=None, mapping_release=None, qc=qc)
            records.append(record)
        provenance = dict(schema_version=1, routing_policy='declared_identifier_species_v1', target_species=target, resources=records,
                          tools={name: importlib.metadata.version(name) for name in ('gseapy', 'pandas', 'numpy')},
                          acquisition_time=None)
        (stage / 'manifest.json').write_text(json.dumps(provenance, indent=2, sort_keys=True) + '\n')
        if output.exists() or output.is_symlink():
            raise FileExistsError(f'Snapshot appeared during staging: {output}')
        os.rename(stage, output)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', required=True, type=Path, help='JSON with explicit requested resources and species evidence')
    parser.add_argument('--output', required=True, type=Path, help='New immutable snapshot directory; existing paths are refused')
    args = parser.parse_args()
    acquire_snapshot(args.manifest, args.output)


if __name__ == '__main__':
    main()
