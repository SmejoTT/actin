from audioop import avg
import json
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# Opening JSON file
f = open('data_sim2.json')
 
# returns JSON object as
# a dictionary
data = json.load(f)
name = data['data_'][0]['first']
X = []
factinY = []
gactinY = []
arpY = []
branchY = []
thymY = []
thymActY = []
for d in data['data_']:
    if d['first'] == "num-agents":
        X = d['second']['x_values']
        factinY = d['second']['y_values']
    if d['first'] == "gactin":
        gactinY = d['second']['y_values']
    if d['first'] == "arp":       
        arpY = d['second']['y_values']
    if d['first'] == "branch":    
        branchY = d['second']['y_values']
    if d['first'] == "thym":
        thymY = d['second']['y_values']
    if d['first'] == "thym_actin":
        thymActY = d['second']['y_values']

X[:] = [x / 1000 for x in X]

def make_figure(xs,ys,name,label1,label2,image_file):
    fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(6,4))
    for i in range(len(xs)):
        if i == 0:
            ax.plot(xs[i], ys[i], color='r',linewidth=2.5,label=label1)
        else:
            ax.plot(xs[i], ys[i], color='b',linewidth=2.5,label=label2)
    ax.set(xlabel='time (min)', ylabel='particle count', title=name)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax.set_xticks(range(0,81,10))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='upper left')
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0,right=80)
    #ax.set_clip_on(False)
    #img = Image.open(image_file)
    #ax[1].axis('off')
    #ax[1] = plt.imshow(img,aspect='auto')
    plt.rcParams['font.size'] = 9
    plt.tight_layout()
    plt.savefig(f'fig_{name}.pdf')
    plt.show()
print(np.average(branchY))
make_figure([X,X],[factinY,gactinY], 'G-actin vs F-actin','F-actin','G-actin','actin.tif')
make_figure([X,X],[branchY,arpY], 'Free Arp vs Bound Arp','Bound Arp2/3','Free Arp2/3','arp.tif')
make_figure([X,X],[thymY,thymActY], 'Thymosin vs Thymosin-Actin','Thymosin','Thymosin-Actin','thym.tif')