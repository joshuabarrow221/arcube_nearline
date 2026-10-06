"""Stream selected segment arrays directly from FLOW, with optional caching.

HTTP byte ranges avoid downloading unrelated hit/packet payloads. No lifetime
attributes or published summaries are consumed. Cached arrays retain exact
values/dtypes and a source URL, ETag, size and per-dataset SHA-256.
"""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import tempfile
from urllib.parse import urlparse

import h5py


DEFAULT_SEGMENTS = ('analysis/beam_rock_muon_segments/data',
                    'analysis/cosmic_muon_segments/data',
                    'analysis/rock_muon_segments/data')


def is_remote(path):
    return str(path).startswith(('https://', 'http://'))


@contextmanager
def open_flow(path):
    """Read local or HTTP FLOW, returning the HDF5 handle and source provenance."""
    if not is_remote(path):
        with h5py.File(path, 'r') as source:
            yield source, dict(input_file=str(Path(path).resolve()), access='local read-only HDF5')
        return
    import aiohttp
    import fsspec
    import requests
    response = requests.head(path, allow_redirects=True, timeout=(10, 30))
    response.raise_for_status()
    size = int(response.headers['Content-Length'])
    etag = response.headers.get('ETag', '')
    if not etag or etag.startswith('W/'):
        raise ValueError('remote FLOW needs a strong ETag for consistent range reads')
    provenance = dict(input_file=path, source_size_bytes=size, source_etag=etag,
                      access='HTTP ranges from FLOW, in-memory block cache')
    filesystem = fsspec.filesystem('http', skip_instance_cache=True,
        client_kwargs={'trust_env': True, 'timeout': aiohttp.ClientTimeout(total=60)})
    with filesystem.open(path, size=size, block_size=262144, cache_type='blockcache',
                         headers={'If-Match': etag, 'Accept-Encoding': 'identity'}) as stream:
        with h5py.File(stream, 'r') as source:
            yield source, provenance
        # fsspec counts requested blocks, including a possibly clipped final block.
        provenance['requested_range_bytes'] = stream.cache.total_requested_bytes


def cache_flow_segments(url, directory, requested_paths=None):
    import aiohttp
    import fsspec
    import requests
    response = requests.head(url, allow_redirects=True, timeout=(10, 30))
    response.raise_for_status()
    size = int(response.headers['Content-Length'])
    etag = response.headers.get('ETag', '')
    if not etag or etag.startswith('W/'):
        raise ValueError('remote FLOW needs a strong ETag for consistent range reads')
    paths = tuple(requested_paths or DEFAULT_SEGMENTS)
    identity = hashlib.sha256((url+'\n'+json.dumps(paths)).encode()).hexdigest()[:16]
    directory = Path(directory)/identity
    directory.mkdir(parents=True, exist_ok=True)
    target = directory/(Path(urlparse(url).path).name+'.segments.h5')
    if target.exists():
        with h5py.File(target, 'r') as cache:
            if (cache.attrs.get('source_flow_url') == url and cache.attrs.get('source_etag') == etag
                    and cache.attrs.get('source_size_bytes') == size and cache.attrs.get('complete', False)):
                return str(target)
    filesystem = fsspec.filesystem('http', skip_instance_cache=True,
        client_kwargs={'trust_env': True, 'timeout': aiohttp.ClientTimeout(total=60)})
    fd, name = tempfile.mkstemp(suffix='.h5.tmp', dir=directory)
    os.close(fd)
    try:
        with filesystem.open(url, size=size, block_size=262144, cache_type='blockcache',
                             headers={'If-Match': etag, 'Accept-Encoding': 'identity'}) as stream:
            with h5py.File(stream, 'r') as source, h5py.File(name, 'w') as cache:
                manifest = []
                for path in paths:
                    if path not in source:
                        if requested_paths:
                            raise KeyError(f'requested segment dataset missing: {path}')
                        continue
                    array = source[path][:]
                    cache.create_dataset(path, data=array)
                    manifest.append(dict(path=path, shape=list(array.shape), dtype=str(array.dtype),
                                         sha256=hashlib.sha256(array.tobytes()).hexdigest()))
                if not manifest:
                    raise KeyError('no selected segments in remote FLOW; full event selection requires a local FLOW file')
                cache.attrs.update(source_flow_url=url, source_etag=etag, source_size_bytes=size,
                    requested_range_bytes=stream.cache.total_requested_bytes, source_datasets=json.dumps(manifest),
                    complete=True, contents='Exact segment arrays from FLOW; no lifetime estimates copied')
        os.chmod(name, 0o644)
        os.replace(name, target)
    finally:
        if os.path.exists(name):
            os.unlink(name)
    return str(target)
