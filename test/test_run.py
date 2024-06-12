import water_shuffle
import pytest
import os
import mdtraj as mdt

rootdir = os.path.dirname(os.path.abspath('__file__'))
trajfile = os.path.join(rootdir, 'test/test.nc')
topfile = os.path.join(rootdir, 'test/test.prmtop')

def test_run():
    t = mdt.load(trajfile, top=topfile)
    std1 = t.xyz.std(axis=0).mean()
    permuted = water_shuffle.run(t, 'HOH', 5)
    std2 = permuted.xyz.std(axis=0).mean()
    assert std2 < std1

