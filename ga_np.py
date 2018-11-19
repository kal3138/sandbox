# -*- coding: utf-8 -*-

import random
import math
import numpy as np

from ga import GA_method


class SGA_method_np(GA_method):

    def __crossover_SPX(self, selected, n_gene):

        n = n_gene + 1
        if len(selected) >= n:
            parents = random.sample(selected, n)
        else:
            parents = [random.choice(selected) for _ in range(n)]
        gene = np.array([p[1] for p in parents])
        gr = np.mean(gene, axis=0)
        eps = np.sqrt(n_gene + 2)
        sr = gr + eps * (gene - gr)
        cr = np.empty([n, n_gene])
        cr[0, :] = sr[0, :]
        rr = np.random.rand
        for i in range(1, n):
            cr[i, :] = math.pow(rr(), 1 / i) * (
                sr[i - 1, :] - sr[i, :] + cr[i - 1, :])

        child = cr.tolist()

        return child
