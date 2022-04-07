import json
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import scipy.io as sio

def make_figure(name,x,ysMy,ysMonk,title,label1,label2,pos):
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(6,4),sharey=True)
    for i in range(len(ysMy)):
        if i == 0:
            ax[0].plot(x, ysMy[i], color='r',linewidth=2.5,label=label1)
        else:
            ax[0].plot(x, ysMy[i], color='b',linewidth=2.5,label=label2)
    for i in range(len(ysMonk)):
        if i == 0:
            ax[1].plot(x, ysMonk[i], color='r',linewidth=2.5,label=label1)
        else:
            ax[1].plot(x, ysMonk[i], color='b',linewidth=2.5,label=label2)
    ax[0].set(xlabel='time (min)', ylabel='particle count', title='My implementation')
    ax[1].set(xlabel='time (min)', title="Monk's implementation")
    
    for i in range(2):
        ax[i].set_xticks(range(0,81,10))
        ax[i].spines['top'].set_visible(False)
        ax[i].spines['right'].set_visible(False)
        #ax[i].axvline(x=43,linestyle='--',c='orange',label="LTP end")
        #ax[i].axvline(x=40,linestyle='--',c='g',label='LTP start')
        ax[i].set_ylim(bottom=0)
        ax[i].set_xlim(left=0,right=80)
    ax[0].legend(loc=pos)
    plt.rcParams['font.size'] = 9
    plt.tight_layout()
    plt.savefig(f'fig_{name}.pdf')
    plt.show()

# Opening JSON file
f = open('data_sim.json')
 
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

BradData = sio.loadmat('ActinData.mat')['ActinData'][0][0]
print(BradData)
make_figure("actin",X,[factinY,gactinY],[np.transpose(BradData[7]),np.transpose(BradData[8])] ,True,'F-actin','G-actin','upper left')
make_figure("arp",X,[branchY,arpY],[np.transpose(BradData[16]),np.transpose(BradData[17])] ,False,'Bound Arp2/3','Free Arp2/3','center right')
make_figure("thymosin",X,[thymY,thymActY],[np.transpose(BradData[21]),np.transpose(BradData[22])] ,True,'Thymosin','Thymosin-Actin','lower center')