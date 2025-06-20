import argparse
import os
import utils
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Run SPADE pipeline")
    parser.add_argument('-s', '--st', type=str, help='Input ST expr file path (required)')
    parser.add_argument('-r', '--sc', type=str, help='Input SC expr file path (required)')
    parser.add_argument('-a', '--sc_anno', type=str, help='Input SC annotation file path (required)')
    parser.add_argument('-c', '--st_coord', type=str, default='', help='Input ST coordinate file path. If not provided, will not enable CCC refinement.')
    parser.add_argument('-o', '--outputpath', type=str, default='./', help='Output file path (default: output.txt)')
    parser.add_argument('--LR', type=str, default='./scriabin_LR_OmniPath.txt', help='The LR database file path (default: ./scriabin_LR_OmniPath.txt)')
    parser.add_argument('--TF', type=str, default='./dorothea.rds', help='The TF database file path (default: ./dorothea.rds)')
    parser.add_argument('--useAllGenes', action='store_true', help='If set, use all genes instead of Marker Genes (default: False)')
    

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    gene = []
    if not args.useAllGenes:
        os.system(f'Rscript codes/FindMarkerGenes.R {args.sc} {args.anno}')
        with open(args.sc[:-4] + '_markers.txt') as g:
            for line in g:
                aa = line.strip()
                gene.append(aa)

    X = pd.read_table(
        args.sc,
        index_col=0,
        sep=",",
    )
    Y = pd.read_table(
        args.st,
        index_col=0,
        sep=",",
    )
    if args.sc_anno:
        os.system(f'Rscript codes/CalcDist.R {args.st_coord}')
        os.system(f'Rscript codes/CalcCCI.R {args.st}')
        anno = pd.read_table(
            args.sc_anno,
            index_col=0,
            sep=",",
        )
        I = pd.read_table(
            args.sc[:-4] + '_CCI.csv',
            index_col=0,
            sep=",",
        )
        D = pd.read_table(
            args.sc[:-4] + '_dist.csv',
            index_col=0,
            sep=",",
        )
        utils.RunAlgorithm(X, Y, anno, I, D, args.outputpath, gene)
    else:
        utils.RunAlgorithm(X, Y, None, None, None, args.outputpath, gene)
    
    os.system(f'Rscript codes/CalcCCI.R {args.outputpath}/Ypredicted.csv')
    os.system(f'Rscript codes/CalcCCCSpots.R {args.outputpath}/Ypredicted_CCI.csv ' + args.sc[:-4] + '_dist.csv {args.outputpath}/CCCScore.csv')
