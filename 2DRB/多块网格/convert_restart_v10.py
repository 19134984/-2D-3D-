"""Convert a connected-ring v10 checkpoint to compact named-array v11 storage.

Usage: python convert_restart_v10.py INPUT.bin OUTPUT.bin
The input, latest.meta and history files are never modified. Choose a new output
path, retain the matching history files, and point latest.meta to the converted file.
Requires NumPy. Binary assumptions match the solver: little-endian int32/float64.
"""
from pathlib import Path
import argparse
import struct
import numpy as np


def convert(source, destination):
    source, destination = Path(source), Path(destination)
    data = np.memmap(source, mode='r', dtype='u1')
    if len(data) < 272 or bytes(data[:16]).rstrip() != b'MB2DRESTART0010':
        raise ValueError('Expected a complete v10 checkpoint header')
    head = np.frombuffer(data, dtype='<i4', count=9, offset=16).copy()
    nx, ny, ratio, left, right, bottom, top, overlap, count = map(int, head)
    if min(nx, ny, ratio) < 1 or count != (1 if ratio == 1 else 2):
        raise ValueError('Unexpected v10 geometry')
    steady = bool(struct.unpack_from('<i', data, 52)[0])
    offset = 272
    coarse_geometry = np.frombuffer(data, dtype='<i4', count=7, offset=offset)
    ni, nj, *_, nh = map(int, coarse_geometry)
    if min(ni, nj) < 1 or nh != (0 if ratio == 1 else 2):
        raise ValueError('Invalid coarse geometry/history')
    coarse_end = offset+84+ni*nj*(20+20*(nh+1)+3*steady)*8
    cuts = []
    fine_arrays = []
    if ratio > 1:
        if coarse_end+84 > len(data):
            raise ValueError('Truncated coarse state')
        geom = np.frombuffer(data, dtype='<i4', count=7, offset=coarse_end)
        coord = np.frombuffer(data, dtype='<f8', count=7, offset=coarse_end+28)
        if list(geom) != [nx, ny, 1, nx, 1, ny, 0] or list(coord) != [0, 0, 1, 0, nx, 0, ny]:
            raise ValueError('Expected full rectangular fine-ring v10 storage')
        nl, nr = left+overlap*ratio, right+overlap*ratio
        nb, nt = bottom+overlap*ratio, top+overlap*ratio
        if min(nx-nl-nr, ny-nb-nt) < 1:
            raise ValueError('Compact layout requires a positive central fine-storage hole')
        cuts = [(0, nl, 0, ny), (nx-nr, nx, 0, ny),
                (nl, nx-nr, 0, nb), (nl, nx-nr, ny-nt, ny)]
        offset = coarse_end+84
        components = [9, 5]+[1]*6+[1]*6+[9, 5]+[1]*(3*steady)
        expected_end = offset+nx*ny*sum(components)*8
        if expected_end != len(data):
            raise ValueError('Truncated state or unexpected trailing checkpoint fields')
        for component_count in components:
            a = np.ndarray((nx, ny, component_count), dtype='<f8', buffer=data,
                           offset=offset, order='F')
            fine_arrays.append(a)
            offset += a.size*8
    elif coarse_end != len(data):
        raise ValueError('Unexpected single-grid checkpoint size')
    head[-1] = 1 if ratio == 1 else 5
    # Exclusive creation prevents replacing either an existing output or the input.
    with destination.open('xb') as output:
        output.write(b'MB2DRESTART0011'.ljust(16))
        output.write(head.tobytes())
        output.write(data[52:coarse_end])
        for x0, x1, y0, y1 in cuts:
            ni, nj = x1-x0, y1-y0
            output.write(struct.pack('<7i', ni, nj, 1, ni, 1, nj, 0))
            output.write(struct.pack('<7d', x0, y0, 1, x0, x1, y0, y1))
            for field in fine_arrays:
                output.write(field[x0:x1, y0:y1, :].tobytes(order='F'))
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    output = convert(args.source, args.destination)
    print(f'Converted: {output}')
    print('Keep the matching NuRe/convergence history; update latest.meta to this filename when ready to resume.')


if __name__ == '__main__':
    main()
