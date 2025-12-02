class LSA:
    def __init__(self, n):
        self.n = n
        self.adj = [[] for _ in range(n)]
        self.E = set()
        self.T = [set() for _ in range(n)]
        self.Tedge = set()
        # self.Tpaths

    def add_edge(self, u, v):
        self.adj[u].append(v)
        self.adj[v].append(u)
        if (u < v):
            self.E.add((u, v))
        else:
            self.E.add((v, u))

    def buildT(self):
        for i in range(self.n):
            self.T[i] = set()
        visited = [False] * self.n
        queue = [0]
        visited[0] = True
        while queue:
            u = queue.pop(0)
            for v in self.adj[u]:
                if not visited[v]:
                    visited[v] = True
                    self.T[u].add(v)
                    self.T[v].add(u)
                    if (u < v):
                        self.Tedge.add((u, v))
                    else:
                        self.Tedge.add((v, u))
                    queue.append(v)

    def buildMST(self):     # if weighted graph
        print("Not implemented yet")

    def lsa(self):
        n = self.n
        self.buildT()

        while (True):   # running phases
            k = 0
            for i in range(n):
                k = max(k, len(self.T[i]))

            while (True):   # running subphases
                mxd = 0
                for i in range(n):
                    mxd = max(mxd, len(self.T[i]))
                if (mxd < k):
                    break

                Dk = set()
                Dk1 = set()
                for i in range(n):
                    if (len(self.T[i]) == k):
                        Dk.add(i)
                    elif (len(self.T[i]) == k-1):
                        Dk1.add(i)

                F = set()
                for (u,v) in self.Tedge:
                    if (u in Dk or u in Dk1 or v in Dk or v in Dk1):
                        # self.Tedge.remove((u,v))
                        F.add((u,v))

                C = dict()
                cid = 0
                for i in range(n):
                    C[i] = -1
                for i in range(n):
                    if (C[i] == -1):
                        C[i] = cid
                        queue = [i]
                        while queue:
                            u = queue.pop(0)
                            for v in self.T[u]:
                                if ((u,v) not in F or (v,u) in F):
                                    if (C[v] == -1):
                                        C[v] = cid
                                        queue.append(v)
                        cid += 1

                Reducible = dict()
                
                while (True):
                    edges_added = []
                    for (v,w) in self.E:
                        if ((C[v] != C[w]) and v not in Dk and v not in Dk1 and w not in Dk and w not in Dk1):
                            edges_added.append((v,w))
                    if (len(edges_added) == 0):
                        return
                    
                    for (v,w) in edges_added:
                        Cycle = []
                        parent = dict()
                        found = False
                        queue = [v]
                        parent[v] = -1
                        while queue and not found:
                            u = queue.pop(0)
                            for x in self.T[u]:
                                if (x not in parent):
                                    parent[x] = u
                                    if (x == w):
                                        found = True
                                        break
                                    queue.append(x)

                        if not found:
                            return
                        
                        cur = w
                        while (cur != -1):
                            Cycle.append(cur)
                            cur = parent[cur]
                        Cycle.reverse()
                        
                        CiDk = []
                        CiDk1 = []
                        for u in Cycle:
                            if u in Dk:
                                CiDk.append(u)
                            elif u in Dk1:
                                CiDk1.append(u)
                        
                        if (len(CiDk) == 0):
                            for u in CiDk1:
                                Reducible[u] = (v,w)
                                Dk1.remove(u)
                            # update F and C
                            for u in CiDk1:
                                for i in range(n):
                                    edge = (u,i)
                                    if i < u:
                                        edge = (i,u)
                                    if edge in F:
                                        F.remove(edge)
                                        for j in range(n):
                                            if C[j] == C[i]:
                                                C[j] = C[u]
                            continue

                        # u = CiDk[0]
                        def local_move(u):
                            if u not in Reducible:
                                return False
                            (v,w) = Reducible[u]
                            local_move(v)
                            local_move(w)

                            Cycle = []
                            parent = dict()
                            found = False
                            queue = [v]
                            parent[v] = -1
                            while queue and not found:
                                x = queue.pop(0)
                                for y in self.T[x]:
                                    if (y not in parent):
                                        parent[y] = x
                                        if (y == u):
                                            found = True
                                            break
                                        queue.append(y)

                            if not found:
                                return False
                            
                            nxt = parent[u]
                            self.Tedge.remove((min(u,nxt), max(u,nxt)))
                            self.T[u].remove(nxt)
                            self.T[nxt].remove(u)
                            self.Tedge.add(min(v,w), max(v,w))
                            self.T[v].add(w)
                            self.T[w].add(v)

                            