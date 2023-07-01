#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import scipy.sparse as sp
import random as rd


# In[2]:


class def_episode:
    def __init__(self, evnts, edge_set):#, lastevnt, pred, succ):
        self.freq = 0
        self.evnts = evnts.copy()
        self.edges = edge_set.copy()
#         self.lastevnt = lastevnt
#         self.pred = pred.copy()
#         self.succ = succ.copy()

class def_NFA:
    def __init__(self):
        self.s = set()
        self.S = list([self.s])
        self.F = -1
        self.D = {}

class def_DFA:
    def __init__(self):
        self.s = {0}
        self.S = list([self.s])
        self.F = list()
        self.D = {}


# In[3]:


def Construct_NFA(alpha):
    alpha_size = len(alpha.evnts)
    pi = {}
    for e in alpha.evnts:
        pi[e] = set()
    for (e1,e2) in alpha.edges:
        pi[e2].add(e1)
    NFA = def_NFA()
    n_states = 1
    for Q in NFA.S:
        NFA.D[str(Q)] = {}
        if len(Q) != alpha_size:
            for ev in alpha.evnts-Q:
                if not (alpha.evnts-Q).intersection(pi[ev]):
                    Q_new = Q.copy()
                    Q_new.add(ev)
                    if Q_new in NFA.S:
                        NFA.D[str(Q)][ev] = {NFA.S.index(Q_new)}
                    else:
                        NFA.S.append(Q_new)
                        NFA.D[str(Q)][ev] = {len(NFA.S)-1}
#                     NFA.D[str(Q)][ev] = [Q_new]
#                     if Q_new not in NFA.S:
#                         NFA.S.append(Q_new)
            NFA.D[str(Q)]['def'] = {NFA.S.index(Q)}
        else:
            NFA.D[str(Q)]['def'] = set()
            NFA.F = NFA.S.index(Q)
    for ev in alpha.evnts:
        if not alpha.evnts.intersection(pi[ev]):
            NFA.D[str(NFA.s)][ev].add(0)
    return NFA


# In[4]:


def Construct_DFA(NFA):
    DFA = def_DFA()
    for state in DFA.S:
        if (NFA.F in state):
            DFA.F.append(state)
        DFA.D[str(state)] = {}
        W = set({})
        for substate in state:
            Q = str(NFA.S[substate])
            W = W.union(set(NFA.D[Q].keys()))
        for ev in W:
            new_state = set()
            for substate in state:
                Q = str(NFA.S[substate])
                if ev in NFA.D[Q].keys():
                    new_state = new_state.union(NFA.D[Q][ev])
#                     DFA.D[str(state)][ev] = DFA.D[str(state)][ev].union(NFA.D[Q][ev])
                else:
                    new_state = new_state.union(NFA.D[Q]['def'])
#                     DFA.D[str(state)][ev] = DFA.D[str(state)][ev].union(NFA.D[Q]['def'])
            DFA.D[str(state)][ev] = new_state.copy()
            if new_state not in DFA.S:
                        DFA.S.append(new_state)
        keys = list(DFA.D[str(state)].keys())
        for k in range(len(keys)):
            if keys[k] != 'def' and DFA.D[str(state)][keys[k]] == DFA.D[str(state)]['def']:
                del DFA.D[str(state)][keys[k]]
            else:
                k += 1
    return DFA


# In[5]:


episodes = list([])
episodes.append(def_episode(set({'1','2','3','4'}),set({('1','2'),('1','3'),('1','4'),('2','3'),('2','4'),('3','4')})))
episodes.append(def_episode(set({'5','6','7','8'}),set({('5','6'),('5','7'),('5','8'),('6','7'),('6','8'),('7','8')})))
NFA = list()
DFA = list()



final_states = list()

for episode in episodes:
    NFA.append(Construct_NFA(episode))
    DFA.append(Construct_DFA(NFA[-1]))
    
    final_states.append(set())
    for k in range(len(DFA[-1].S)):
        if (NFA[-1].F in DFA[-1].S[k]):
            final_states[-1].add(k)


# In[6]:


M = 50
T = 10000

qn = 0.5
qe = 0

pe = 1/M

curr_state = list([])
for k in range(len(episodes)):
    curr_state.append({0})

noise_flag = True
count = [0]*len(episodes)

datastream = list([])
for t in range(1,T+1):
    num = rd.uniform(0,1)
    
    if noise_flag:
        if num > qn:
            noise_flag = False
    else:
        if num > qe:
            noise_flag = True
    
    if noise_flag:
        datastream.append((str(rd.sample(range(1,M+1),1)[0]),t))
    else:
        ep_id = rd.sample(range(len(episodes)),1)[0]
        state = curr_state[ep_id]
        nextevents = list(DFA[ep_id].D[str(state)].keys())
        if NFA[ep_id].F not in state:
            nextevents.remove('def')
            new_event = rd.sample(nextevents,1)[0]
            datastream.append((new_event,t))
            curr_state[ep_id] = DFA[ep_id].D[str(state)][new_event]
        else:
            count[ep_id] += 1
            probs = list([])
            for nxtevnt in nextevents:
                if nxtevnt == 'def':
                    probs.append(1-pe*(len(nextevents)-1))
                else:
                    probs.append(pe)
            new_event = rd.choices(nextevents, probs, k=1)[0]
            if new_event != 'def':
                datastream.append((new_event,t))
                curr_state[ep_id] = DFA[ep_id].D[str(state)][new_event]
            else:
                datastream.append((str(rd.sample(set(range(1,M+1))-episodes[ep_id].evnts,1)[0]),t))
                curr_state[ep_id] = DFA[ep_id].D[str(state)][new_event]
print('No of occurremces of episodes embedded:', count)


# In[7]:


name = 'Data/Datastream_M_'+ str(M) + '_n_'+ str(T) + '_N_4_emb_' + str(len(episodes)) + '_q1_' + str(qn) + '_q2_' + str(qe) + '_set_4.txt'
with open(name,'w') as f:

    for (e,t) in datastream:
        f.write(e + ',' + str(t) + '\n')
    


# In[ ]:




