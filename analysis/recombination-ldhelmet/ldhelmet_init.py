'''
ldhelmet_init.py - run LDhelmet on aligned regions
generated by align_mt_fasta_maf.py

run from project root dir
'''

import argparse
import os
import glob
import subprocess
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(description = 'ldhelmet_init.py - initial LDhelmet run',
                                     usage = 'ldhelmet_int.py [options]')

    parser.add_argument('-d', '--directory', required = True,
                        type = str, help = 'Directory containing files.')
    parser.add_argument('-l', '--ldhelmet', required = True,
                        type = str, help = 'LDhelmet executable path')
    parser.add_argument('-o', '--outdir', required = True,
                        type = str, help = 'Folder for all final outfiles')

    args = parser.parse_args()

    return [args.directory, args.ldhelmet, args.outdir]

def run_ldhelmet(alignment_file, ldhelmet_path, directory, outdir):

    if not os.path.exists(directory + 'temp'):
        os.mkdir(directory + 'temp')
    if not directory.endswith('/'):
        directory += '/'
    if not outdir.endswith('/'):
        outdir += '/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    coordinates = alignment_file.split('/')[-1].split('.')[0]

    find_confs = '{ldhelmet} find_confs \
    --num_threads 10 \
    --window_size 50 \
    --output_file {directory}temp/{coordinates}.conf \
    {alignment_file}'.format(
        ldhelmet = ldhelmet_path, 
        directory = directory,
        coordinates = coordinates,
        alignment_file = alignment_file)

    table_gen = '{ldhelmet} table_gen \
    --num_threads 10 \
    --conf_file {directory}temp/{coordinates}.conf \
    --theta 0.01 \
    --rhos 0.0 0.1 10.0 1.0 100.0 \
    --output_file {directory}temp/{coordinates}.lk'.format(
        ldhelmet = ldhelmet_path,
        coordinates = coordinates,
        directory = directory)

    pade = '{ldhelmet} pade \
    --num_threads 10 \
    --conf_file {directory}temp/{coordinates}.conf \
    --theta 0.01 \
    --output_file {directory}temp/{coordinates}.pade'.format(
        ldhelmet = ldhelmet_path,
        coordinates = coordinates,
        directory = directory)

    rjmcmc = '{ldhelmet} rjmcmc \
    --num_threads 10 \
    --window_size 50 \
    --seq_file {alignment_file} \
    --lk_file {directory}temp/{coordinates}.lk \
    --pade_file {directory}temp/{coordinates}.pade \
    --num_iter 1000000 --burn_in 100000 \
    --block_penalty 50 \
    --output_file {directory}temp/{coordinates}.post'.format(
        ldhelmet = ldhelmet_path,
        directory = directory,
        coordinates = coordinates,
        alignment_file = alignment_file)

    post_to_text = '{ldhelmet} post_to_text \
    --mean \
    --perc 0.025 \
    --perc 0.50 \
    --perc 0.975 \
    --output_file {outdir}{coordinates}.txt \
    {directory}temp/{coordinates}.post'.format(
        ldhelmet = ldhelmet_path,
        outdir = outdir,
        coordinates = coordinates,
        directory = directory)

    cmds = [find_confs, table_gen, pade, rjmcmc, post_to_text]
    print('LDhelmet run for {f}...'.format(f = coordinates))
    for cmd in cmds:
        child = subprocess.Popen(cmd, stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE, shell = True)
        stdout, stderr = child.communicate()
        cmd_name = cmd.split(' ')[1]
        print('Completed {cmd}.'.format(cmd = cmd_name))

    # cleanup
    print('Done.')
    print('Outfile written to {outfile}'.format(outfile = outdir + coordinates + '.txt'))
    print('Removing temp files...')
    temp_dir = '{directory}temp/*'.format(directory = directory)
    for fname in glob.glob(temp_dir):
        os.remove(fname)

def main():
    directory, ldhelmet, outdir = args()

    # create filename list to iterate over
    if not directory.endswith('/'):
        directory += '/'
    fnames = glob.glob(directory + '*fasta')
    
    for alignment in tqdm(fnames):
        run_ldhelmet(alignment, ldhelmet, directory, outdir)

if __name__ == '__main__':
    main()

        
