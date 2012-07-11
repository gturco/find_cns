import subprocess


def create_env(dir_name):
    """creates virtenv with pip and python 2.7"""
    create_env = subprocess.call(['virtualenv','--distribute', dir_name,'--python=python2.7'])
    #assert: create_env == 0

def pip_install(dir_name):
    """pip install packages to virenv dir"""
    subprocess.call(['pip', 'install', '-E', dir_name, 'numpy'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'processing'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'shapely'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'pyfasta'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'scipy'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'Cython'])
    subprocess.call(['pip', 'install', '-E', dir_name, 'biopython'])



def git_install(dir_name):
    """downloads git scripts to virenv bin"""
    #flatfeature = subprocess.call(['pip', 'install', '-E', dir_name, 'git+https://github.com/brentp/flatfeature.git'])
    #if flatfeature != 0:
    #    print "download git to contiune"
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
