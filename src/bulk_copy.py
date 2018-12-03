import os
from shutil import copyfile
from glob import glob

def main():

    savdir = '../results/planet_yield_plots/all_ideas_one_directory/'
    names = ['idea_1_SNE','idea_2_SNSNS','idea_3_SNNSN','idea_4_EC3PO']
    outdirs = [os.path.join('../results/planet_yield_plots/',dirname)
               for dirname in names]

    for outdir in outdirs:
        fpaths = glob(os.path.join(outdir,'*png'))

        if len(fpaths)==0:
            print('got no files in {}'.format(outdir))

        else:
            for fpath in fpaths:
                newname = fpath.split('/')[-1].rstrip('.png')
                newname += '_'+fpath.split('/')[-2]
                newname += '.png'
                outpath = os.path.join(savdir,newname)
                copyfile(fpath, outpath)
                print('{} --> {}'.format(fpath, outpath))

if __name__=="__main__":

    main()
