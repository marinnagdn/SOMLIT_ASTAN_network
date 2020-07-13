import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import sys
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from colour import Color
import mplcursors
import itertools
from collections import defaultdict
import seaborn as sns
from math import floor,ceil
import pygraphviz
import cmocean
import igraph as ig
import random
import matplotlib.patches as mpatches
import re
import ast
from collections import Counter
import matplotlib.gridspec as gd
from venn import venn


def flat_list(l):
	# go to list of lists to a simple flat list #
	try:
		return(sum(l, []))
	except TypeError:
		return(sum(l, ()))

#####################


def normalize(df):
    result = df.copy()
    for feature_name in df.columns:
        max_value = df[feature_name].max()
        min_value = df[feature_name].min()
        result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
    return result
    

######################

def rgb_to_hex(rgb):
    # for network display #
    return '%02x%02x%02x' % rgb

####################

def mycolormap(index=None):
    colors = ["#6D93EE","#E01D22","#FABD28","#96AF00",
    "#804DFF","#FF7749","#47240F","#4CCDA4","#B60085",
    "#3B9A3C","#2B5985","#FF429B","#5FFFF8","#20FC86",
    "#9600A5","#B1E7FF","#80FF00","#6B4C81","#989898",
    "#0F7C6A","#BB4C71",'#B4FF9B',"#EAD788","#FFC6F6",
    "#FF00C2","#E2A4FF","#FF4A00","#800037","#1F45FF"]*5
    if index == None:
        return colors
    else:
        return(colors[index])