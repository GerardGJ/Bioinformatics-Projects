#!/usr/bin/env python

import sys
import re


class Seq:
    def __init__(self):
        self.id = ""
        self.seq = ""
        self.features = ""

    def __str__(self):
        n = 60
        seq_chunks = [self.seq[i : i + n] for i in range(0, len(self.seq), n)]
        seq_chunks = "\n".join(seq_chunks)
        return f">{self.id}\n{seq_chunks}"


class Node:
    def __init__(self):
        self.name = ""
        self.distance = 0.0
        self.bootstrap = 0.0
        self.left = None
        self.right = None
        self.parent = None


def read_seq(fname):
    seq_list = {}
    records = []
    record_id = -1
    with open(fname) as f:
        for line in f:
            if line.startswith(">"):
                records.append(Seq())
                match_obj = re.match(r">(\S*)\s*(.*)", line)
                record_id += 1
                if match_obj:
                    records[record_id].id = match_obj.group(1)
                    records[record_id].features = match_obj.group(2)
                records[record_id].seq = ""
            else:
                records[record_id].seq += line.strip()

    for i in range(record_id + 1):
        records[i].seq = re.sub(r"[ \n\t\r]", "", records[i].seq)
        seq_list[records[i].id] = records[i].seq
    return seq_list


def read_matrix(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    matrix = {}  # will be a dictionary of dictionaries matrix[aa1][aa2]

    for line in lines[:-1]:
        splitted = line.split()
        aa1 = splitted[0]
        matrix[aa1] = {}
        aa_alphabet = list(matrix.keys())

        for aa2, numb in zip(aa_alphabet, splitted[1:]):
            matrix[aa1][aa2] = int(numb)
            matrix[aa2][aa1] = int(numb)
    return matrix


def align(seq_in, Igroup, Jgroup, matrix, gep):
    seq = seq_in.copy()
    lenI = len(seq[Igroup[0]])
    lenJ = len(seq[Jgroup[0]])

    score_mat = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]
    tb = [[0 for x in range(lenJ + 1)] for y in range(lenI + 1)]

    for i in range(0, lenI + 1):
        score_mat[i][0] = i * gep
        tb[i][0] = 1

    for j in range(0, lenJ + 1):
        score_mat[0][j] = j * gep
        tb[0][j] = -1

    for i in range(1, lenI + 1):
        for j in range(1, lenJ + 1):
            score = 0
            nsubst = 0
            for ni in range(0, len(Igroup)):
                for nj in range(0, len(Jgroup)):
                    a1 = seq[Igroup[ni]][i - 1]
                    a2 = seq[Jgroup[nj]][j - 1]
                    if a1 != "-" and a2 != "-":
                        score += int(matrix[a1.upper()][a2.upper()])
                        nsubst += 1
            if nsubst > 0:
                score /= nsubst

            Sub = score_mat[i - 1][j - 1] + score
            Del = score_mat[i][j - 1] + gep
            Ins = score_mat[i - 1][j] + gep

            if Sub > Del and Sub > Ins:
                score_mat[i][j] = Sub
                tb[i][j] = 0
            elif Del > Ins:
                score_mat[i][j] = Del
                tb[i][j] = -1
            else:
                score_mat[i][j] = Ins
                tb[i][j] = 1

    i = lenI
    j = lenJ
    new_seq = {}
    for ni in range(0, len(Igroup)):
        new_seq[Igroup[ni]] = ""
    for nj in range(0, len(Jgroup)):
        new_seq[Jgroup[nj]] = ""

    while (i == 0 and j == 0) != 1:
        if tb[i][j] == 0:
            i -= 1
            j -= 1
            for ni in range(0, len(Igroup)):
                new_seq[Igroup[ni]] += seq[Igroup[ni]][i]
            for nj in range(0, len(Jgroup)):
                new_seq[Jgroup[nj]] += seq[Jgroup[nj]][j]

        elif tb[i][j] == -1:
            j -= 1
            for ni in range(0, len(Igroup)):
                new_seq[Igroup[ni]] += "-"
            for nj in range(0, len(Jgroup)):
                new_seq[Jgroup[nj]] += seq[Jgroup[nj]][j]

        elif tb[i][j] == 1:
            i -= 1
            for ni in range(0, len(Igroup)):
                new_seq[Igroup[ni]] += seq[Igroup[ni]][i]
            for nj in range(0, len(Jgroup)):
                new_seq[Jgroup[nj]] += "-"

    for ni in range(0, len(Igroup)):
        seq[Igroup[ni]] = new_seq[Igroup[ni]][::-1]
    for nj in range(0, len(Jgroup)):
        seq[Jgroup[nj]] = new_seq[Jgroup[nj]][::-1]
    return seq, int(score_mat[lenI][lenJ])

def join(tree,name,i,j):
    tx = tree.pop(max(i,j))
    ty = tree.pop(min(i,j))
    nx = name.pop(max(i,j))
    ny = name.pop(min(i,j))
    name.append(nx+ny)
    tree.append([tx,ty])
    return tree, name

def tree_creator(tree_list):
    nodes = []
    a = Node()
    a.left = len(nodes)+1
    nodes.append(a)
    nodes = tree_creator_aux(tree_list,0,nodes)
    return nodes


def tree_creator_aux(tree_list,parent,nodes):
    for x in range(len(tree_list)):
        a = Node()
        a.parent = parent
        num_node = len(nodes)
        if len(tree_list) == 2 and x == 1:
            nodes[a.parent].right = num_node
        if type(tree_list[x][0]) == str:
            a.name = tree_list[x][0]
            nodes.append(a)
        else:
            a.left = len(nodes)+1
            nodes.append(a)
            nodes = tree_creator_aux(tree_list[x],num_node,nodes)
    return nodes

#Read the file with the fasta format sequences and store the in a diccionary
seq_dic = read_seq(sys.argv[1])
matrix = read_matrix(sys.argv[2])
seq_names = seq_dic.keys()
gep = int(sys.argv[3])
name_list = []
tree_list = []

for record in seq_dic:
    name_list.append([record])
    tree_list.append([record])


while len(name_list) > 1:
    good_align, maximum = align(seq_dic,name_list[0],name_list[1],matrix,gep)
    maxi = 0
    maxj = 1
    
    for i in range(len(name_list)):
        for j in range(i+1,len(name_list)):
            new_align, num = align(seq_dic,name_list[i],name_list[j],matrix,gep)
            if num > maximum:
                maximum = num
                good_align = new_align
                maxi = i
                maxj = j
    seq_dic = good_align
    print('Joining nodes', tree_list[maxi], 'and', tree_list[maxj])
    tree_list, name_list = join(tree_list,name_list,maxi,maxj)

nodes = tree_creator(tree_list)

for name in seq_dic:
    print('Sequence:', name)
    print(seq_dic[name])
for x in range(len(nodes)):
    print('Node',x,'=',nodes[x].right,nodes[x].left,nodes[x].name)
