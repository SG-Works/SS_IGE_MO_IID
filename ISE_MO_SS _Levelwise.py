#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy.stats as st
import math
import sympy as smp
import random as rd
import matplotlib.pyplot as plt


# In[2]:


class def_event:
    def __init__(self,E,t):
        self.evnt = E
        self.time = t

class def_serial_episode:
    def __init__(self, evnts_ordr):
        self.evnt = evnts_ordr[:]
#         self.length = len(evnts_ordr)
        self.freq = 0

class def_window:
    def __init__(self, st_time, end_time):
        self.ts = st_time
        self.te = end_time


# In[3]:


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


# In[4]:


def MO_Span_Recursive(ep_p,ep_p_MO_list,Tx,FT):
    global F1
    F = list()
    closed = 1
    if len(ep_p.evnt) <= 3:
        for E in F1:
            if E not in ep_p.evnt:
                evnts_temp = ep_p.evnt.copy()
                evnts_temp.append(E)
                ep = def_serial_episode(evnts_temp)
                ep_MO_list = Compute_MO_Span_List(ep_p_MO_list,E)
                ep.freq = len(ep_MO_list)
                if ep.freq == ep_p.freq:
                    closed = 0
                if ep.freq >= FT:
                    F += MO_Span_Recursive(ep,ep_MO_list,Tx,FT)
#     if closed == 1:
    F.append(ep_p)
    return(F)


# In[5]:


def Compute_MO_Span_List(ep_p_MO_list,E):
    p_size = len(ep_p_MO_list)
    E_MO_list = MO_list[evnt_id[E]]
    E_size = len(E_MO_list)
    ep_MO_list = list()
    p1=0
    p2=0
    while (p1<p_size and p2<E_size):
        while p2 < E_size and E_MO_list[p2].ts <= ep_p_MO_list[p1].te :
            p2 += 1
        if p2 < E_size and E_MO_list[p2].te-ep_p_MO_list[p1].ts < Tx:
            while(p1 < p_size and ep_p_MO_list[p1].te < E_MO_list[p2].ts):
                p1 += 1
            ep_MO_list.append(def_window(ep_p_MO_list[p1-1].ts,E_MO_list[p2].te))
        else:
            p1 += 1
    return(ep_MO_list)


# In[6]:


def Find_Closed(Freq_ep):
    hasher = list()
    hashtable = list()
    for ep in Freq_ep:
        fn = ep.freq
        flg = 1
        for k in range(len(hasher)):
            if fn == hasher[k]:
                flg = 0
                hashtable[k].append(ep)
                break
        if flg == 1:
            hasher.append(fn)
            hashtable.append(list([ep]))

    Freq_closed_ep = list()
    for k in range(len(hasher)):
        len_hash_k = len(hashtable[k])
        p = 0
        while p < len_hash_k:
            q = 0
            if q == p:
                q += 1
            while q < len_hash_k:
                len_p = len(hashtable[k][p].evnt)
                len_q = len(hashtable[k][q].evnt)
                if len_p < len_q:
                    sub_ep = 0
                    pi = 0
                    qi = 0
                    while(1):
                        if hashtable[k][p].evnt[pi] == hashtable[k][q].evnt[qi]:
                            pi += 1
                            if pi == len_p:
                                sub_ep = 1
                                break
                        qi += 1
                        if qi == len_q:
                            break
                    if sub_ep == 1:
                        hashtable[k].remove(hashtable[k][p])
                        len_hash_k = len(hashtable[k])
                        q -= 1
                        break
                q += 1
                if q == p:
                    q += 1
            if q >= len_hash_k:
                p += 1

    for k in range(len(hasher)):
        for l in range(len(hashtable[k])):
            Freq_closed_ep.append(hashtable[k][l])

    return Freq_closed_ep


# In[7]:


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


# In[8]:


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


# In[9]:


def GetTransitionMatrix(DFA,prob):
    T0 = list()
    for state in DFA.S:
        T0.append([0]*len(DFA.S))
        for ev in DFA.D[str(state)].keys():
            if ev == 'def':
                continue
            T0[-1][DFA.S.index(DFA.D[str(state)][ev])] += prob[ev]
        T0[-1][DFA.S.index(DFA.D[str(state)]['def'])] = 1 - sum(T0[-1])
    return(T0)


# In[10]:


def GetStats(T0,final_states):
    
    u = smp.symbols('u')
    Dphi = smp.eye(len(T0))
    for final in final_states:
        Dphi[final,final] = u

    Tu = smp.ImmutableSparseMatrix(T0)*smp.ImmutableSparseMatrix(Dphi)
#     smp.pprint(T0)
#     smp.pprint(Dphi)
#     smp.pprint(Tu)

    lamda = smp.symbols('lamda')
#     Qn = (smp.ImmutableMatrix(lamda*smp.eye(len(T0)))-Tu).det()
    
    Qn = Tu.charpoly(lamda)
#     print(Qn)
#     print(Qn.as_expr())

    dQdu = smp.diff(Qn.as_expr(),u).subs([(u,1),(lamda,1)])
    dQdl = smp.diff(Qn.as_expr(),lamda).subs([(u,1),(lamda,1)])
    d2Qdu2 = smp.diff(smp.diff(Qn.as_expr(),u),u).subs([(u,1),(lamda,1)])
    d2Qdl2 = smp.diff(smp.diff(Qn.as_expr(),lamda),lamda).subs([(u,1),(lamda,1)])
    d2Qdudl = smp.diff(smp.diff(Qn.as_expr(),u),lamda).subs([(u,1),(lamda,1)])

    lambda1 = -dQdu/dQdl
    lambda11 = -(d2Qdu2 + 2*lambda1*d2Qdudl + lambda1**2*d2Qdl2)/dQdl
    # print(lambda1,lambda11)

    MX1 = lambda1
    VX1 = lambda11 + lambda1 - lambda1**2

#     print(MX1,VX1)
#     MX1 = 0
#     VX1 = 1
    return(MX1, VX1)


# In[11]:


# ISE Length2

evs = list(['1','2'])
ISE2 = def_episode(set(evs),set({('1','2')}))

NFA = Construct_NFA(ISE2)
DFA = Construct_DFA(NFA)

final_states = set()
for k in range(len(DFA.S)):
    if (NFA.F in DFA.S[k]):
        final_states.add(k)

Prob_syms2 = {}
probs2 = smp.symbols('pe:'+str(len(evs)))
for k in range(len(evs)):
    Prob_syms2[evs[k]] = probs2[k]


T0 = GetTransitionMatrix(DFA,Prob_syms2)
MX1_2, VX1_2 = GetStats(T0,final_states)


# In[12]:


# ISE Length3

evs = list(['1','2','3'])
ISE3 = def_episode(set(evs),set({('1','2'),('1','3'),('2','3')}))

NFA = Construct_NFA(ISE3)
DFA = Construct_DFA(NFA)

final_states = set()
for k in range(len(DFA.S)):
    if (NFA.F in DFA.S[k]):
        final_states.add(k)

Prob_syms3 = {}
probs3 = smp.symbols('pe:'+str(len(evs)))
for k in range(len(evs)):
    Prob_syms3[evs[k]] = probs3[k]

# Prob_values = {}
# ppp = [0.035, 0.035, 0.035]
# for k in range(len(evs)):
#     Prob_values[evs[k]] = ppp[k]

T0 = GetTransitionMatrix(DFA,Prob_syms3)
MX1_3, VX1_3 = GetStats(T0,final_states)



# mmm = MX1_3.subs([(probs[k],Prob_values[evs[k]]) for k in range(len(evs))])
# vvv = VX1_3.subs([(probs[k],Prob_values[evs[k]]) for k in range(len(evs))])

# T0 = GetTransitionMatrix(DFA,Prob_values)
# MX1_3, VX1_3 = GetStats(T0,final_states)


# In[13]:


# ISE Length4

evs = list(['1','2','3','4'])
ISE4 = def_episode(set(evs),set({('1','2'),('1','3'),('1','4'),('2','3'),('2','4'),('3','4')}))

NFA = Construct_NFA(ISE4)
DFA = Construct_DFA(NFA)

final_states = set()
for k in range(len(DFA.S)):
    if (NFA.F in DFA.S[k]):
        final_states.add(k)

Prob_syms4 = {}
probs4 = smp.symbols('pe:'+str(len(evs)))
for k in range(len(evs)):
    Prob_syms4[evs[k]] = probs4[k]

# Prob_values = {}
# ppp = [0.035, 0.035, 0.035, 0.035]
# for k in range(len(evs)):
#     Prob_values[evs[k]] = ppp[k]

# T0 = GetTransitionMatrix(DFA,Prob_values)
# for row in T0:
#     print(row)
# MX1_4, VX1_4 = GetStats(T0,final_states)

T0 = GetTransitionMatrix(DFA,Prob_syms4)
# for row in T0:
#     smp.pprint(row)
MX1_4, VX1_4 = GetStats(T0,final_states)

# mmm = MX1_4.subs([(probs[k],Prob_values[evs[k]]) for k in range(len(evs))])
# vvv = VX1_4.subs([(probs[k],Prob_values[evs[k]]) for k in range(len(evs))])


# In[98]:


evnt_strm = list()
# print('Enter the link of the txt file containing the event stream (for example: ./Data/Datastream_M_50_n_10000_N_4_emb_2_q1_1_q2_0.txt)')
# name = input()
name = './Data/Datastream_M_50_n_10000_N_4_emb_2_q1_1_q2_0.txt'
with open(name,'r') as f:
    for line in f:
        entry = line.split(',')
        evnt_strm.append(def_event(entry[0],int(entry[1][:-1])))

M = len(evnt_strm)
alph = list()
for m in range(len(evnt_strm)):
    A = evnt_strm[m].evnt
    flg = 1
    for a in alph:
        if a == A:
            flg = 0
            break
    if flg == 1:
        alph.append(A)
LA = len(alph) 

A = alph[0]
evnt_id = {A: 0}
count = 1
for A in alph[1:]:
    evnt_id[A] = count
    count += 1

MO_list = list()
for i in range(LA):
    MO_list.append(list())

Freq_ep = list()
F1 = alph
F1_ep = list()
for A in alph:
    F1_ep.append(def_serial_episode(list([A])))



Tx = 75
FT = 70
# print('Enter the Expiry time "Tx":')
# Tx = float(input())
# print('Enter the frequency threshold "FT":')
# FT = float(input())

for event in evnt_strm:
    E = event.evnt
    t = event.time
    F1_ep[evnt_id[E]].freq += 1
    MO_list[evnt_id[E]].append(def_window(t,t))


for E in F1:
    if F1_ep[evnt_id[E]].freq >= FT:
        Freq_ep += MO_Span_Recursive(F1_ep[evnt_id[E]],MO_list[evnt_id[E]],Tx,FT)

# print('The following are the frequent episodes discovered')
# for ep in Freq_ep:
#     print(ep.evnt ,':', ep.freq)
    
Freq_closed_ep = Find_Closed(Freq_ep)

print('Total closed frequent episodes: ', len(Freq_closed_ep))
count_FE_2 = 0
count_FE_3 = 0
count_FE_4 = 0
for k in range(len(Freq_closed_ep)):
    if len(Freq_closed_ep[k].evnt) == 2:
        count_FE_2 += 1
    if len(Freq_closed_ep[k].evnt) == 3:
        count_FE_3 += 1
    if len(Freq_closed_ep[k].evnt) == 4:
        count_FE_4 += 1
print(count_FE_2)
print(count_FE_3)
print(count_FE_4)


# In[93]:


serial_episodes = list()
for ep in Freq_closed_ep:
    if len(ep.evnt) >=2 and len(ep.evnt) <=4:
        serial_episodes.append(ep)

n = 10000
Prob_Ep_Events = {}
for E in F1:
    Prob_Ep_Events[E] = len(MO_list[evnt_id[E]])/n

mean_estimate = list()
var_estimate = list()

for episode in serial_episodes:
    if len(episode.evnt) == 2:
        mmm = MX1_2.subs([(probs2[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        vvv = VX1_2.subs([(probs2[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        mean_estimate.append(mmm)
        var_estimate.append(vvv)
    if len(episode.evnt) == 3:
        mmm = MX1_3.subs([(probs3[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        vvv = VX1_3.subs([(probs3[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        mean_estimate.append(mmm)
        var_estimate.append(vvv)
    if len(episode.evnt) == 4:
        mmm = MX1_4.subs([(probs4[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        vvv = VX1_4.subs([(probs4[k],Prob_Ep_Events[episode.evnt[k]]) for k in range(len(episode.evnt))])
        mean_estimate.append(mmm)
        var_estimate.append(vvv)
    print(mmm,vvv)


# In[96]:


SSTh = list()
c = 0.0001
for k in range(len(serial_episodes)):
    value = mean_estimate[k]*n + math.sqrt(var_estimate[k]*n)*st.norm.ppf(1-c)
    SSTh.append(value)

count = 0
for k in range(len(serial_episodes)):
    print(serial_episodes[k].evnt, ' : ', serial_episodes[k].freq, ' , ', SSTh[k])
    if serial_episodes[k].freq > SSTh[k]:
        count += 1


# In[97]:


print(len(serial_episodes),count)
count_FE_2 = 0
count_FE_3 = 0
count_FE_4 = 0
count_SS_2 = 0
count_SS_3 = 0
count_SS_4 = 0
for k in range(len(serial_episodes)):
    if len(serial_episodes[k].evnt) == 2:
        count_FE_2 += 1
        if serial_episodes[k].freq > SSTh[k]:
            count_SS_2 += 1
    if len(serial_episodes[k].evnt) == 3:
        count_FE_3 += 1
        if serial_episodes[k].freq > SSTh[k]:
            count_SS_3 += 1
    if len(serial_episodes[k].evnt) == 4:
        count_FE_4 += 1
        if serial_episodes[k].freq > SSTh[k]:
#             print(serial_episodes[k].evnt, ' : ', serial_episodes[k].freq, ' , ', SSTh[k])
            count_SS_4 += 1
print(count_SS_2)
print(count_SS_3)
print(count_SS_4)

