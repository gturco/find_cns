import subprocess
import sys

def create_env(dir_name):
    """creates virtenv with pip and python 2.7"""
    print >>sys.stderr, 'Make sure all requirements of INSTALL file are downloaded before running!!!'
    create_env = subprocess.call(['virtualenv','--distribute', dir_name,'--python=python2.7'])
    #assert: create_env == 0

def pip_install(dir_name):
    """pip install packages to virenv dir"""
    numpy = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'numpy'])
    if numpy != 0: raise Exception('can not download numpy READ REQUIREMENTS and TRY again!')
    processing = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'processing'])
    if processing != 0: raise Exception('can not download processing READ REQUIREMENTS and TRY again!')
    shapely = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'shapely'])
    if shapely != 0: raise Exception('can not download shapely READ REQUIREMENTS and TRY again!')
    pyfasta = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'pyfasta'])
    if pyfasta != 0: raise Exception('can not download pyfasta READ REQUIREMENTS and TRY again!')
    scipy = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'scipy'])
    if scipy != 0: raise Exception('can not download scipy READ REQUIREMENTS and TRY again!')
    Cython = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'Cython'])
    if Cython != 0: raise Exception('can not download Cython READ REQUIREMENTS and TRY again!')
    biopython = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'biopython'])
    if biopython != 0: raise Exception('can not download biopython READ REQUIREMENTS and TRY again!')
    

def git_install(dir_name):
    """downloads git scripts to virenv bin"""
    print >>sys.stderr, 'Be patient, takes a long time to download'
    flatfeature = subprocess.call(['{0}/bin/pip'.format(dir_name), 'install', 'git+https://github.com/brentp/flatfeature.git'])
    if flatfeature != 0:
        raise Exception("Download git to contiune")
    quota = subprocess.Popen(['git', 'clone','https://github.com/tanghaibao/quota-alignment.git'],cwd=r'{0}/bin/'.format(dir_name)).wait()
    bcbb = subprocess.Popen(['git', 'clone', 'https://github.com/chapmanb/bcbb.git'],cwd=r'{0}/bin/'.format(dir_name)).wait()


def setup_install(dir_name):
    """installs setup install files to virenv directory"""
    subprocess.Popen(['../../python2.7','setup.py','install'],cwd=r'{0}/bin/bcbb/gff/'.format(dir_name))
    subprocess.Popen(['../../../../{0}/bin/python2.7'.format(dir_name),'setup.py','install'],cwd=r'pipeline/coann/brents_bpbio/biostuff/').wait()
    subprocess.Popen(['../../../../../{0}/bin/python2.7'.format(dir_name),'setup.py','install'],cwd=r'pipeline/coann/brents_bpbio/blasttools/blast_misc/').wait()
    subprocess.Popen(['../../../../../{0}/bin/python2.7'.format(dir_name),'setup.py','install'],cwd=r'pipeline/coann/brents_bpbio/scripts/bblast/').wait()
    co_anno = subprocess.Popen(['../../../../{0}/bin/python2.7'.format(dir_name),'setup.py','install'],cwd=r'pipeline/coann/brents_bpbio/co-anno/').wait()





create_env('cns_pipeline')
pip_install('cns_pipeline')
git_install('cns_pipeline')
setup_install('cns_pipeline')
