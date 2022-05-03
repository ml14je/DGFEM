#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joseph Elmes: NERC-Funded PhD Researcher in Applied Mathematics
University of Leeds : Leeds LS2 9JT : ml14je@leeds.ac.uk

Python 3.7
Created on Thu Jan 14 19:26:59 2021
"""
import numpy as np

def sparse_matlab(i, j, v, m, n):
    from scipy.sparse import csr_matrix as sp
    A = sp((m, n))

    try:
        ni = i.shape[1]

    except IndexError:
        ni = 1
        i, j = i[:, None], j[:, None]

    try:
        v.shape[0]

    except AttributeError:
        v = v * np.ones((ni, i.shape[0]))

    try:
        vn2 = v.shape[0]

    except IndexError:
        v = v[:, None]
        vn2 = 1

    assert ni == vn2

    for ii in range(ni):
        B = sp((v[ii, :], (i[:, ii], j[:, ii])), shape=(m, n))
        A += B

    return A

def refine(geo, p, e, t, it, mode='regular'):
    """
        REFINEMESH Refine a triangular mesh.
        [P1,E1,T1]=REFINEMESH(P,E,T) returns a refined version
        of the triangular mesh specified by point matrix P,
        edge matrix E, and triangle matrix T.

        The triangular mesh is given by the mesh data P, E, and T.
        Details can be found under INITMESH.

        [P1,E1,T1,U1]=REFINEMESH(P,E,T,U) refines the mesh and also
        extends the function U to the new mesh by linear interpolation.
        The number of rows in U should correspond to the number of
        columns in P, and U1 will have as many rows as there are
        points in P1. Each column of U is interpolated separately.

        An extra input argument is interpreted as a list of
        subdomains to refine, if it is a row vector, or a list of
        triangles to refine, if it is a column vector.

        The default refinement method is regular refinement, where
        all of the specified triangles are divided into four triangles
        of the same shape. Longest edge refinement, where
        the longest edge of each specified triangle is bisected, can
        be demanded by giving 'longest' as a final parameter. Using
        'regular' as the final parameter results in regular refinement.
        Some triangles outside of the specified set may also be refined,
        in order to preserve the triangulation and its quality.

        See also INITMESH, PDEGEOM

        Note: This function does not support quadratic triangular elements


        A. Nordmark 4-26-94, AN 6-21-94.
        Copyright 1994-2017 The MathWorks, Inc.
    """
    Np, Ne, Nt = p.shape[0], e.shape[0], t.shape[0]
    # For 32 GB RAM, max number of elements is  8 Bytes * 32 * 1024^3/64 Bytes = 4 * 1024^3 elements
    maxsize = 4 * (1024**3)
    indexproblem = Np**2 > maxsize
    print(indexproblem)

    # Find longest side of each triangle
    ls = 3 * np.ones((Nt, 1))
    d1 = (p[t[:, 0], 0] - p[t[:, 0], 1])**2 + (p[t[:, 0], 1]-p[t[:, 1], 1])**2
    d = (p[t[:, 1], 0] - p[t[:, 2], 0])**2 + (p[t[:, 1], 1]-p[t[:, 2], 1])**2
    ls[d > d1] = 1
    d1 = np.maximum(d, d1)
    d = (p[t[:, 2], 0] - p[t[:, 0], 0])**2 + (p[t[:, 2], 1]-p[t[:, 0], 1])**2
    ls[d > d1] = 2
    # Permute so longest side is 3
    ii = ls[:, 0] == 1
    d = t[ii, 0]
    t[ii, 0] = t[ii, 1]
    t[ii, 1] = t[ii, 2]
    t[ii, 2] = d
    ii = ls[:, 0] == 2
    d = t[ii, 0]
    t[ii, 0] = t[ii, 2]
    t[ii, 2] = t[ii, 1]
    t[ii, 1] = d
    itt1 = np.ones((1, Nt))
    itt1[it] = 0
    it = np.where(itt1[0] == 0)[0] # Triangles whos longest side is to be bisected
    it1 = np.where(itt1[0] != 0)[0] # Triangles not yet to be refined

    # Make a connectivity matrix, with edges to be refined.
    # -1 means no point is yet allocated
    ip1 = t[it, 0]
    ip2 = t[it, 1]

    if mode == 'regular':
        ip3 = t[it, 2]

    A = sparse_matlab(ip1[:, None], ip2[:, None], -1, Np, Np)
    if mode == 'regular':
        A += sparse_matlab(ip2[:, None], ip3[:, None], -1, Np, Np)
        A += sparse_matlab(ip3[:, None], ip1[:, None], -1, Np, Np)

    A = -1*((A + A.T)<0)
    newpoints=1
    assert A.shape[0] == Np

    # loop until no additional hanging nodes are introduced
    while newpoints:
        newpoints=0
        n = len(it1)
        ip1, ip2, ip3 = t.T[:, it1]
        assert ip1.shape[0] == n

        m1, m2, m3 = np.zeros(n), np.zeros(n), np.zeros(n)
        # print(it1.shape, m3[it1].shape, ip1[it1].shape, ip2[it1].shape, A.shape)
        m3, m1, m2 = A[ip1, ip2], A[ip2, ip3], A[ip3, ip1]
        ii = np.nonzero(m3)[-1]
        if len(ii) > 0:
            itt1[0, it1[ii]] = 0

        ii = np.nonzero(((m1 | m2) & (~m3)))[-1]
        if len(ii) > 0:
            A += sparse_matlab(ip1[ii], ip2[ii], -1, Np, Np)
            A = -1*((A + A.T)<0)
            newpoints = 1
            itt1[0, it1[ii]] = 0

        it1 = np.nonzero(itt1)[-1] #Triangles not yet fully refined
        it = np.where(itt1 == 0)[-1] #Triangles fully refined

        if not indexproblem:
            ie = A[e[:, 0], e[:, 1]] == -1
        else:
            ie = l_extract(A, e[:, 0], e[:, 1]) == -1

        ie1 = ie == 0 #Edges not to be refined
        ie = ie != 0 #Edges to be refined

        print(ie)
        print(ie1)

        from PDE import pdeigeon
        x, y = pdeigeon(geo, e[4, ie], (e[2, ie]+e[3, ie])/2)
        p1 = np.array([p, np.array([x, y])])


# % Find edges to be refined
# if ~indexproblem
#   ie=full(A(e(1,:)+(e(2,:)-1)*np))==-1;
# else
#   ie=l_extract(A,e(1,:),e(2,:))==-1;
# end

# ie1=find(ie==0);                        % Edges not to be refined
# ie=find(ie);                            % Edges to be refined

# % Get the edge "midpoint" coordinates
# [x,y]=pdeigeom(g,e(5,ie),(e(3,ie)+e(4,ie))/2);
# % Create new points
# p1=[p [x;y]];
# if intp
#   u1=[u;(u(e(1,ie),:)+u(e(2,ie),:))/2];
# end
# ip=(np+1):(np+length(ie));
# np1=np+length(ie);
# % Create new edges
# e1=[e(:,ie1) ...
#         [e(1,ie);ip;e(3,ie);(e(3,ie)+e(4,ie))/2;e(5:7,ie)] ...
#         [ip;e(2,ie);(e(3,ie)+e(4,ie))/2;e(4,ie);e(5:7,ie)]];
# % Fill in the new points
# if ~indexproblem
#   A(e(1,ie)+np*(e(2,ie)-1))=ip;
#   A(e(2,ie)+np*(e(1,ie)-1))=ip;
# else
#   A=l_assign(A,[e(1,ie) e(2,ie)],[e(2,ie) e(1,ie)],[ip ip]);
# end

# % Generate points on interior edges
# [i1,i2]=find(A==-1 & A.'==-1);
# i=find(i2>i1);
# i1=i1(i);
# i2=i2(i);
# p1=[p1 ((p(1:2,i1)+p(1:2,i2))/2)];
# if intp
#   u1=[u1;(u(i1,:)+u(i2,:))/2];
# end
# ip=(np1+1):(np1+length(i));
# np1=np1+length(i);
# % Fill in the new points
# if ~indexproblem
#   A(i1+np*(i2-1))=ip;
#   A(i2+np*(i1-1))=ip;
# else
#   A=l_assign(A,[i1 i2],[i2 i1],[ip ip]);
# end

# % Lastly form the triangles
# ip1=t(1,it);
# ip2=t(2,it);
# ip3=t(3,it);
# if ~indexproblem
#   mp1=full(A(ip2+np*(ip3-1)));
#   mp2=full(A(ip3+np*(ip1-1)));
#   mp3=full(A(ip1+np*(ip2-1)));
# else
#   mp1=l_extract(A,ip2,ip3);
#   mp2=l_extract(A,ip3,ip1);
#   mp3=l_extract(A,ip1,ip2);
# end

# % Find out which sides are refined
# bm=1*(mp1>0)+2*(mp2>0);
# % The number of new triangles
# nt1=length(it1)+length(it)+sum(mp1>0)+sum(mp2>0)+sum(mp3>0);
# t1=zeros(4,nt1);
# t1(:,1:length(it1))=t(:,it1);           % The unrefined triangles
# nnt1=length(it1);
# if isempty(bm)
# 	i = bm;
# else
# 	i=find(bm==3);                          % All sides are refined
# end
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp1(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[mp1(i);mp2(i);mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# if isempty(bm)
# 	i = bm;
# else
# 	i=find(bm==2);                          % Sides 2 and 3 are refined
# end
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(1,it(i));mp3(i);mp2(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp2(i);mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# if isempty(bm)
# 	i = bm;
# else
# 	i=find(bm==1);                          % Sides 3 and 1 are refined
# end
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));mp1(i);mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));mp3(i);mp1(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# if isempty(bm)
# 	i = bm;
# else
# 	i=find(bm==0);                          % Side 3 is refined
# end
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(3,it(i));t(1,it(i));mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);
# t1(:,(nnt1+1):(nnt1+length(i)))=[t(2,it(i));t(3,it(i));mp3(i);t(4,it(i))];
# nnt1=nnt1+length(i);

def l_extract(A,i,j):
    i1, j1 = i.flatten(), j.flatten()
    if len(i1) != len(j1):
        raise TypeError('Incorrect number of elements in i and j.')

    k = np.zeros(len(i1))

    for l in range(len(i1)):
        k[l] = A[i[l], j[l]]

    return k

def l_assign(A, i, j, k):
    i1, j1, k1 = i.flatten(), j.flatten(), k.flatten()
    il = len(i1)
    if il != len(j1) or il != len(k1):
        raise TypeError('Incorrect number of elements in i, j and k.')

    for l in range(il):
        A[i1[l], j1[l]] = k1[l]

    return A



if __name__ == '__main__':
    pass