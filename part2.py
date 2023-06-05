import pandas as pd
import argparse
import seaborn as sns
import os
import glob
import matplotlib.pyplot as plt

def plot_gibbs(input_directory):
    '''
    Function to read scores files and plot the scores distributions
    '''
    all_files = glob.glob(os.path.join(input_directory, "*.csv"))
    # print()
    frames = []
    
    for filename in all_files:
        df = pd.read_csv(filename)
        frames.append(df)
    
    scores = pd.concat(frames, axis=1)
    scores["distance"] = range(0, 21)
    
    return scores
    
def plot_distributions(input_directory, output_directory):
    '''
    Function to plot the scores distributions
    '''
    scores = plot_gibbs(input_directory)
    
    for column in scores.columns:
        if column != "distance":
            plt.figure()
            sns.lineplot(x=scores["distance"], y=scores[column], color='b')
            plt.xlabel("Distance (Å)")
            plt.ylabel("Score")
            plt.title("Pair {}".format(column))
            plt.savefig(os.path.join(output_directory, column))
            plt.show()
            plt.close()


