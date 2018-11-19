# -*- coding: utf-8 -*-

import random
import math
import statistics
import copy
from operator import add, sub
from itertools import chain


class GA(object):

    def __init__(self, eval_func):

        self.EVALUATE_FUNCTION = eval_func  # 評価関数

        try:
            super().__init__()
        except AttributeError:
            pass


class GA_score(object):

    def get_score(self, population):

        score_list = [None] * len(population)
        for i, p in enumerate(population):
            score_list[i] = p['score']

        try:
            super().get_score(population)
        except AttributeError:
            pass

        return score_list

    def calc_score(self, x, eval_func):

        try:
            super().calc_score(x)
        except AttributeError:
            pass

        return eval_func(x)

    def evaluate(self, population, eval_func):

        for p in population:
            if p['score'] is None:
                p['score'] = self.calc_score(p['gene'], eval_func)

        try:
            super().evaluate(population, eval_func)
        except AttributeError:
            pass

        return fitness


class SGA(GA):

    default_value = {
        'crossover_rate': 0.8,
        'mutate_rate': 0.05,
        'elite_rate': 0.2,
        'n_tournament': 4,
        'preserve_elite': True,
        'lethal_gene': False,
        'death_function': None,
        'cross_type': '2point',
        'select_type': 'roulette',
        'mutation_type': 'substitution',
        'alpha': 0.3,
        'generation_model': 'MGG',
        'large_mutation_rate': 0,
        'large_mutation_condition': None
    }

    def __init__(self, n_gene, n_po, n_generation, eval_func):

        super().__init__(eval_func)
        self.N_GENE = n_gene  # 遺伝子長
        self.N_POP = n_po  # 1世代の個体数
        self.N_GENERATION = n_generation  # 最大計算世代数

        # set_paramで変更可能
        # 交配率
        self.CROSSOVER_RATE = SGA.default_value['crossover_rate']
        # 突然変異率
        self.MUTATE_RATE = SGA.default_value['mutate_rate']
        # 全個体数に対するエリート個体の割合
        self.ELITE_RATE = SGA.default_value['elite_rate']
        # トーナメントサイズ
        self.N_TOURNAMENT = SGA.default_value['n_tournament']
        # エリート保存実行フラグ
        self.PRESERVE_ELITE = SGA.default_value['preserve_elite']
        # 致死遺伝子の有無
        self.LETHAL_GENE = SGA.default_value['lethal_gene']
        # 致死関数
        self.DEATH_FUNCTION = SGA.default_value['death_function']
        # 交配方法
        self.CROSS_TYPE = SGA.default_value['cross_type']
        # 選択方法
        self.SELECT_TYPE = SGA.default_value['select_type']
        # 突然変異方法
        self.MUTATION_TYPE = SGA.default_value['mutation_type']
        #
        self.ALPHA = SGA.default_value['alpha']
        # 世代交代モデル
        self.GENERATION_MODEL = \
            SGA.default_value['generation_model']
        # 大変異率
        self.LARGE_MUTATION_RATE = \
            SGA.default_value['large_mutation_rate']
        # 大変異発生条件
        self.LARGE_MUTATION_CONDITION = SGA.default_value[
            'large_mutation_condition']


class SGA_select(object):

    def __select_elite(self, population, n_pop, elite_rate):

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)
        elites = population[:n_pop * elite_rate]

        return elites

    def __select_roulette(self, population, n_pop):

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)
        scores = super().get_score(population)
        scores -= min(scores)
        total = sum(scores)

        roulette = [None] * n_pop
        for i in range(n_pop):
            target = random.uniform(0, total)
            score = 0
            for j, sc in enumerate(scores):
                score += sc
                if score > target:
                    roulette[i] = j
                    break
            else:
                roulette[i] = j
        sr = set(roulette)
        selected = [None] * len(sr)
        j = 0
        for i in sr:
            selected[j] = population[i]
            j += 1

        return selected

    def __select_rank(self, population, n_pop):

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)
        scores = iter([1 / (i + 1) for i in range(n_pop)])
        total = sum(scores)
        rank = [None] * n_pop
        for i in range(n_pop):
            target = random.uniform(0, total)
            score = 0
            for j, sc in enumerate(scores):
                score += sc
                if score > target:
                    rank[i] = j
                    break
            else:
                rank[i] = j
        sr = set(rank)
        selected = [None] * len(sr)
        j = 0
        for i in sr:
            selected[j] = population[i]
            j += 1

        return selected

    def __select_tournament(
            self, population, n_pop, n_tournament):

        n = math.ceil(n_pop / n_tournament)
        selected = [None] * n
        random.shuffle(population)
        for i in range(n):
            t_round = population[
                i * n_tournament:(i + 1) * n_tournament]
            t_round.sort(
                key=lambda x: (x['score'] is not None, x),
                reverse=True)
            selected[i] = t_round[0]
        return selected

    def select_method(
            self,
            population,
            select_type,
            preserve_elite,
            elite_rate=None,
            n_tournament=None,
    ):

        elites = []
        n_pop = len(population)
        if preserve_elite:
            # Get elites
            if elite_rate is None:
                print('elite_rate is None')
                return
            elites = self.__select_elite(
                population, n_pop, elite_rate)

        if select_type == 'elite':
            selected = copy.deepcopy(elites)

        elif select_type == 'roulette':
            selected = self.__select_roulette(population, n_pop)

        elif select_type == 'rank':
            selected = self.__select_rank(population, n_pop)

        elif select_type == 'tournament':
            if n_tournament is None:
                print('n_tournament is None')
                return
            selected = self.__select_tournament(
                population, n_pop, n_tournament)

        else:
            print('Unknown select_type')
            return

        return elites, selected


class SGA_crossover(object):

    def __crossover_2point(self, parent1, parent2):

        length = len(parent1)
        r1 = random.randint(0, length - 1)
        r2 = random.randint(r1 + 1, length)

        child = copy.deepcopy(parent1)
        child[r1:r2] = parent2[r1:r2]

        return child

    def __crossover_uniform(self, parent1, parent2):

        rd = random.random
        child = [
            g2 if rd() < 0.5 else g1
            for g1, g2 in zip(parent1, parent2)]

        return child

    def __crossover_BLX_alpha(self, parent1, parent2, a):

        adr = (
            a * abs(g1 - g2) for g1, g2 in zip(parent1, parent2))
        min_r = (min(g1, g2) for g1, g2 in zip(parent1, parent2))
        max_r = (max(g1, g2) for g1, g2 in zip(parent1, parent2))

        min_cr = map(sub, min_r, adr)
        max_cr = map(add, max_r, adr)

        ru = random.uniform
        child = [ru(r1, r2) for r1, r2 in zip(min_cr, max_cr)]

        return child

    def crossover(self,
                  parent1,
                  parent2,
                  cross_type,
                  gene_type=None,
                  alpha=None):

        if cross_type == '2point':
            child = self.__crossover_2point(parent1, parent2)
        elif cross_type == 'uniform':
            child = self.__crossover_uniform(parent1, parent2)
        elif cross_type == 'BLX' and gene_type == 'real':
            if alpha is None:
                print('alpha is None')
                return
            child = self.__crossover_BLX_alpha(
                parent1, parent2, alpha)
        else:
            print('Unknown crossover_type')
            return

        return child


class SGA_mutate(object):

    def __mutate_substitution(self, parent, gene_type, gene):

        r = random.randint(0, len(parent) - 1)
        child = copy.deepcopy(parent)
        if gene_type == 'binary':
            c = (parent[r] + 1) % 2
        elif gene_type == 'real':
            g_min, g_max = gene
            c = random.uniform(g_min, g_max)
        else:
            print('Unknown gene_type')
            return
        child[r] = c

        return child

    def __mutate_inversion(self, parent):

        length = len(parent)
        r1 = random.randint(0, length - 1)
        r2 = random.randint(r1 + 1, length)
        child = copy.deepcopy(parent)
        child[r1:r2] = child[r2:r1:-1]

        return child

    def __mutate_scramble(self, parent):

        length = len(parent)
        r1 = random.randint(0, length - 1)
        r2 = random.randint(r1 + 1, length)
        child = copy.deepcopy(parent)
        random.shuffle(child[r1:r2])

        return child

    def __mutate_translocation(self, parent):

        length = len(parent)
        r1 = random.randint(0, length - 1)
        r2 = random.randint(r1 + 1, length)
        dr = r2 - r1
        r_trans1 = random.randint(0, length - dr)
        r_trans2 = r_trans1 + dr
        child = copy.deepcopy(parent)
        dr = abs(r1 - r_trans1)
        if r1 > r_trans1:
            child[r_trans1:r_trans2], \
                child[r_trans2:r_trans2 + dr] \
                = child[r1:r2], child[r_trans1:r1]
        elif r1 < r_trans1:
            child[r_trans1:r_trans2], \
                child[r_trans1:r_trans1 + dr] \
                = child[r1:r2], child[r2:r_trans2]

        return child

    def __mutate_perturbation(self, parent, d_gene):

        dg_min, dg_max = d_gene
        r = int(random.uniform(0, len(parent)))
        dg = random.uniform(dg_min, dg_max)
        child = copy.deepcopy(parent)
        child[r] += dg

        return child

    def mutate(self,
               parent,
               mutation_type,
               gene_type='binary',
               gene=None,
               d_gene=None):

        if mutation_type == 'substitution':
            if gene_type == 'real' and gene is None:
                print('gene is None')
                return
            child = self.__mutate_substitution(
                parent, gene_type, gene)

        elif mutation_type == 'inversion':
            child = self.__mutate_inversion(parent)

        elif mutation_type == 'scramble':
            child = self.__mutate_scramble(parent)

        elif mutation_type == 'translocation':
            child = self.__mutate_translocation(parent)

        if mutation_type == 'perturbation' \
                and gene_type == 'real':
            child = self.__mutate_perturbation(parent, d_gene)

        else:
            print('Unknown mutation_type')
            return

        return child


class SGA_generation(
        GA_score, SGA_select, SGA_crossover, SGA_mutate):

    def __gc_normal(
            self, n_pop, gene, d_gene, func, elites, selected,
            lethal_gene, crossover_rate, mutate_rate, alpha,
            type_param):

        population = elites[:]
        next_pop = []
        cross_type, mutation_type, gene_type = type_param
        eval_func, death_func = func
        die = False

        rr = random.random
        rs = random.sample
        # crossover
        while len(population) < n_pop:
            if rr() < crossover_rate:
                parent1, parent2 = rs(selected, 2)
                child = super().crossover(
                    parent1['gene'], parent2['gene'], cross_type,
                    gene_type, alpha)
                if lethal_gene:
                    die = death_func(child)
                if not die:
                    population.append(
                        {'score': None, 'gene': child})

        # mutate
        for p in population:
            if rr() < mutate_rate:
                mutant = super().mutate(
                    p['gene'], mutation_type, gene_type, gene,
                    d_gene)
                if lethal_gene:
                    die = death_func(mutant)
                    while die:
                        mutant = super().mutate(
                            p['gene'], mutation_type, gene_type,
                            gene, d_gene)
                        die = death_func(mutant)
                next_pop.append({'score': None, 'gene': mutant})
            else:
                next_pop.append(p)

        # Evaluate individual
        super().evaluate(next_pop, eval_func)

        return next_pop

    def __mgg(
            self, population, n_pop, gene, d_gene, func,
            lethal_gene, crossover_rate, mutate_rate, alpha,
            gene_type):

        rr = random.random
        rc = random.choice
        next_pop = population
        eval_func, death_func = func
        if rr() < crossover_rate:
            p1, p2 = random.sample(range(n_pop), 2)
            parent1, parent2 = population[p1], population[p2]
            gene_p1, gene_p2 = parent1['gene'], parent2['gene']
            cross_types = ['2point', 'uniform']
            if gene_type == 'real':
                cross_types.append('BLX')

            family = []
            die = False

            # crossover
            while len(family) < n_pop:
                cross_type = rc(cross_types)
                child = super().crossover(
                    gene_p1, gene_p2, cross_type, gene_type,
                    alpha)
                if lethal_gene:
                    die = death_func(child)
                if not die:
                    family.append({'score': None, 'gene': child})
            family.extend([parent1, parent2])

            # Evaluate individual
            super().evaluate(family, eval_func)

            # Select individual
            # Select elite
            family.sort(
                key=lambda x: (x['score'] is not None, x),
                reverse=True)
            individual1 = family[0]
            # Select roulette
            selected = super().select_method(
                family[1:], 'roulette', False)

            individual2 = random.choice(selected)
            next_pop[p1], next_pop[p2] = individual1, individual2

        # mutate
        mutation_types = [
            'substitution', 'inversion', 'scramble',
            'translocation'
        ]
        if gene_type == 'real':
            mutation_types.append('perturbation')
        for r, p in enumerate(next_pop):
            if rr() < mutate_rate:
                mutation_type = rc(mutation_types)
                if mutation_type != 'BLX':
                    mutant = super().mutate(
                        p['gene'], mutation_type, gene_type,
                        gene, d_gene)
                    if lethal_gene:
                        die = death_func(mutant)
                        while die:
                            mutant = super().mutate(
                                p['gene'], mutation_type,
                                gene_type, gene, d_gene)
                            die = death_func(mutant)
                    next_pop[r] = {'score': None, 'gene': mutant}

        # Evaluate individual
        super().evaluate(next_pop, eval_func)

        return next_pop

    def generation_change(
            self, population, n_pop, generation_model, eval_func,
            crossover_rate, mutate_rate, gene_type, gene, d_gene,
            preserve_elite, elite_rate, n_tornament, select_type,
            cross_type, mutation_type, lethal_gene, death_func,
            alpha):

        func = eval_func, death_func
        type_param = cross_type, mutation_type, gene_type
        if generation_model == 'simple':
            elites, selected = super().select_method(
                population, select_type,
                preserve_elite,
                elite_rate, n_tournament)
            # Cross and mutate
            population = self.__gc_normal(
                n_pop, gene, d_gene, func, elites, selected,
                lethal_gene, crossover_rate, mutate_rate, alpha,
                type_param)
        elif generation_model == 'MGG':
            population = self.__mgg(
                population, n_pop, gene, d_gene, func,
                lethal_gene, crossover_rate, mutate_rate,
                alpha, gene_type)

        return population


class SGA_method(SGA_generation):

    def __get_population(self, n_pop, n_gene, gene, gene_type):

        if gene_type == 'binary':
            for _ in range(n_pop):
                population = [
                    random.randint(0, 1) for _ in range(n_gene)]
                yield population
        elif gene_type == 'real':
            g_min, g_max = gene
            for _ in range(n_pop):
                population = [
                    random.uniform(g_min, g_max)
                    for _ in range(n_gene)
                ]
                yield population
        else:
            print('Unknown gene_type')
            return

    def GA_main(
            self, n_pop, n_gene, n_generation, eval_func,
            generation_model, crossover_rate, mutate_rate,
            gene_type, gene, d_gene, select_type, cross_type,
            mutation_type, pe, elite_rate, n_tornament,
            lethal_gene, death_func, large_mutation,
            large_mutation_condition, alpha):

        population = [
            {'score': None, 'gene': p}
            for p in self.__get_population(
                n_pop, n_gene, gene, gene_type)
        ]
        super().evaluate(population, eval_func)

        for g in range(n_generation):
            print('Generation: ' + str(g + 1))

            if large_mutation > 0:
                f_lmr = large_mutation_condition()
                if f_lmr:
                    mutate_rate = large_mutation

            # Cross and mutate
            population = super().generation_change(
                population, n_pop, generation_model, eval_func,
                crossover_rate, mutate_rate, gene_type, gene,
                d_gene, pe, elite_rate, n_tornament, select_type,
                cross_type, mutation_type, lethal_gene,
                death_func, alpha)

            score = super().get_score(population)

            print('min: ', min(score))  # 最も低い評価を表示
            print('max: ', max(score))  # 最も高い評価を表示
            print('ave: ', statistics.mean(score))  # 評価の平均値を表示

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)

        return population


class GA_single(SGA):

    def __init__(
            self, n_gene, n_po, n_generation, eval_func,
            gene_type, g_min, g_max, dg_min, dg_max):

        super().__init__(n_gene, n_po, n_generation, eval_func)
        if gene_type == 'real':
            self.GENE_MIN = g_min  # 遺伝子の最小値
            self.GENE_MAX = g_max  # 遺伝子の最大値
            self.d_gene_MIN = dg_min  # 遺伝子の微小変化最小値
            self.d_gene_MAX = dg_max  # 遺伝子の微小変化最大値

    def set_param(self,
                  cr=None,
                  mu=None,
                  el=None,
                  nt=None,
                  pe=None,
                  lg=None,
                  df=None,
                  ct=None,
                  st=None,
                  mt=None,
                  al=None,
                  gm=None,
                  lmr=None,
                  lmc=None):

        if cr is not None:
            self.CROSSOVER_RATE = cr
        if mu is not None:
            self.MUTATE_RATE = mu
        if el is not None:
            self.ELITE_RATE = el
        if nt is not None:
            self.N_TOURNAMENT = nt
        if pe is not None:
            self.PRESERVE_ELITE = pe
        if lg is not None:
            self.LETHAL_GENE = lg
        if df is not None:
            self.DEATH_FUNCTION = df
        if ct is not None:
            self.CROSS_TYPE = ct
        if st is not None:
            self.SELECT_TYPE = st
        if mt is not None:
            self.MUTATION_TYPE = mt
        if al is not None:
            self.ALPHA = al
        if gm is not None:
            self.GENERATION_MODEL = gm
        if lmr is not None:
            self.LARGE_MUTATION_RATE = lmr
        if lmc is not None:
            self.LARGE_MUTATION_CONDITION = lmc


class GA_binary(GA_single, SGA_method):

    def __init__(self, n_gene, n_po, n_generation, eval_func=sum):

        super().__init__(
            n_gene, n_po, n_generation, eval_func, 'binary', _, _,
            _, _)

    def main(self):

        population = super().GA_main(
            self.N_POP, self.N_GENE, self.N_GENERATION,
            self.EVALUATE_FUNCTION, self.GENERATION_MODEL,
            self.CROSSOVER_RATE, self.MUTATE_RATE, 'binary', _, _,
            self.SELECT_TYPE, self.CROSS_TYPE, self.MUTATION_TYPE,
            self.PRESERVE_ELITE, self.ELITE_RATE,
            self.N_TOURNAMENT, self.LETHAL_GENE,
            self.DEATH_FUNCTION, self.LARGE_MUTATION_RATE,
            self.LARGE_MUTATION_CONDITION, _)
        print('best')
        print(population[0])  # 最も高い評価の個体を表示


class GA_real(GA_single, SGA_method):

    def __init__(
            self, n_gene, n_po, n_generation, eval_func, g_min,
            g_max, dg_min, dg_max):

        super().__init__(
            n_gene, n_po, n_generation, eval_func, 'real', g_min,
            g_max, dg_min, dg_max)

    def main(self):

        gene = (self.GENE_MIN, self.GENE_MAX)
        d_gene = (self.d_gene_MIN, self.d_gene_MAX)

        population = super().GA_main(
            self.N_POP, self.N_GENE, self.N_GENERATION,
            self.EVALUATE_FUNCTION, self.GENERATION_MODEL,
            self.CROSSOVER_RATE, self.MUTATE_RATE, 'real', gene,
            d_gene, self.SELECT_TYPE, self.CROSS_TYPE,
            self.MUTATION_TYPE, self.PRESERVE_ELITE,
            self.ELITE_RATE, self.N_TOURNAMENT, self.LETHAL_GENE,
            self.DEATH_FUNCTION, self.LARGE_MUTATION_RATE,
            self.LARGE_MUTATION_CONDITION, self.ALPHA)
        print('best')
        print(population[0])  # 最も高い評価の個体を表示


class GA_island(SGA):

    def __init__(
            self, n_gene, n_po, n_generation, n_ge_i, n_island,
            mig_pr, eval_func, g_min, g_max, dg_min, dg_max,
            gene_type):

        if gene_type == 'binary':
            self.ISLAND = [
                GA_binary(n_gene, n_po, n_ge_i, eval_func)
                for i in range(n_island)
            ]
        elif gene_type == 'real':
            self.ISLAND = [
                GA_real(
                    n_gene, n_po, n_ge_i, eval_func, g_min, g_max,
                    dg_min, dg_max) for i in range(n_island)
            ]
        self.N_ISLAND = n_island  # 島数(並列分散数)
        self.N_GENERATION_MAX = n_generation  # 最大計算世代数
        self.MIGRATION_RATE = mig_pr  # 移住割合

    def set_param(self,
                  num,
                  cr=None,
                  mu=None,
                  el=None,
                  nt=None,
                  pe=None,
                  lg=None,
                  df=None,
                  ct=None,
                  st=None,
                  mt=None,
                  al=None,
                  gm=None):

        if cr is None:
            cr = (super().default_value['crossover_rate']
                  for _ in range(self.N_ISLAND))
        if mu is None:
            mu = (super().default_value['mutate_rate']
                  for _ in range(self.N_ISLAND))
        if el is None:
            el = (super().default_value['elite_rate']
                  for _ in range(self.N_ISLAND))
        if nt is None:
            nt = (super().default_value['n_tournament']
                  for _ in range(self.N_ISLAND))
        if pe is None:
            pe = (super().default_value['preserve_elite']
                  for _ in range(self.N_ISLAND))
        if lg is None:
            lg = (super().default_value['lethal_gene']
                  for _ in range(self.N_ISLAND))
        if df is None:
            df = (super().default_value['death_function']
                  for _ in range(self.N_ISLAND))
        if ct is None:
            ct = (super().default_value['cross_type']
                  for _ in range(self.N_ISLAND))
        if st is None:
            st = (super().default_value['select_type']
                  for _ in range(self.N_ISLAND))
        if mt is None:
            mt = (super().default_value['mutation_type']
                  for _ in range(self.N_ISLAND))
        if al is None:
            al = (super().default_value['alpha']
                  for _ in range(self.N_ISLAND))
        if gm is None:
            gm = (super().default_value['generation_model']
                  for _ in range(self.N_ISLAND))

        if isinstance(num, list):
            for i in num:
                j = i - 1
                self.ISLAND[j].set_param(
                    cr[j], mu[j], el[j], nt[j], pe[j], lg[j],
                    ct[j], st[j], mt[j], al[j])
        else:
            n = num - 1
            self.ISLAND[n].set_param(
                cr[n], mu[n], el[n], nt[n], pe[n], lg[n], ct[n],
                st[n], mt[n], al[n])

    def main(self):

        # 各島でGA処理
        for g in range(self.N_GENERATION_MAX):
            print('Generation all: ' + str(g + 1))

            population = [
                self.ISLAND[i].main()
                for i in range(self.N_ISLAND)
            ]

            # 各島の優秀な個体を移住
            for il in range(self.N_ISLAND):
                n_migration = int(
                    len(population[il]) * self.MIGRATION_RATE)
                migrant = population[il][:n_migration]
                del population[il][:n_migration]
                im = il
                while im == il:
                    im = random.randint(0, self.N_ISLAND - 1)
                population[im].extend(migrant)

        population = list(chain.from_iterable(population))
        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)
        print('best')
        print(population[0])  # 最も高い評価の個体を表示

        try:
            super().main()
        except AttributeError:
            pass


class GA_island_binary(GA_island):

    def __init__(self,
                 n_gene,
                 n_po,
                 n_generation,
                 n_ge_i,
                 mig_pr,
                 n_island,
                 eval_func=sum):

        super().__init__(
            n_gene, n_po, n_generation, n_ge_i, n_island, mig_pr,
            eval_func, _, _, _, _, 'binary')


class GA_island_real(GA_island):

    def __init__(
            self, n_gene, n_po, n_generation, n_ge_i, mig_pr,
            n_island, eval_func, g_min, g_max, dg_min, dg_max):

        super().__init__(
            n_gene, n_po, n_generation, n_ge_i, n_island, mig_pr,
            eval_func, g_min, g_max, dg_min, dg_max, 'real')


class PfGA(GA):

    def __init__(
            self, n_gene, n_generation, eval_func, lethal_gene,
            death_func):

        super().__init__(eval_func)
        self.N_GENE = n_gene  # 遺伝子長
        self.N_GENERATION = n_generation  # 最大計算世代数
        self.LETHAL_GENE = lethal_gene  # 致死遺伝子の有無
        self.DEATH_FUNCTION = death_func  # 致死関数


class PfGA_method(GA_score):

    def __get_population(
            self, n_pop, n_gene, g_min, g_max, gene_type):

        if gene_type == 'binary':
            ri = random.randint
            population = [
                [ri(0, 1) for _ in range(n_gene)]
                for _ in range(n_pop)]
        elif gene_type == 'real':
            ru = random.uniform
            population = [
                [ru(g_min, g_max) for _ in range(n_gene)]
                for _ in range(n_pop)]
        if n_pop == 1:
            population = list(chain.from_iterable(population))

        return population

    def __select_parents(self, population):

        return random.sample(population, 2)

    def __crossover(
            self, parent1, parent2, lethal_gene, death_func):

        length = len(parent1)
        child = copy.deepcopy(parent1)

        ri = random.randint
        rs = random.sample
        for c in range(2):
            die = True
            while die:
                n = ri(1, length)
                p = iter(rs(list(range(length)), n))
                for r in p:
                    child[r] = parent2[r]

                if lethal_gene:
                    die = death_func(child)
                else:
                    die = False

            if c == 1:
                child1 = {'score': None, 'gene': child}
            else:
                child2 = {'score': None, 'gene': child}

        return child1, child2

    def __mutate(
            self, child, lethal_gene, death_func, g_min, g_max,
            gene_type):

        length = len(child)
        die = True
        while die:
            mutant = copy.deepcopy(child)
            n = random.randint(0, length)
            if n > 0:
                r = iter(random.sample(list(range(length)), n))
                if gene_type == 'binary':
                    for i in r:
                        mutant[i] = (child[i] + 1) % 2
                elif gene_type == 'real':
                    for i in r:
                        mutant[i] = random.uniform(g_min, g_max)
            if lethal_gene:
                die = death_func(child)
            else:
                die = False

        return {'score': None, 'gene': mutant}

    def generation_change(
            self, population, eval_func, n_gene, lethal_gene,
            death_func, g_min, g_max, gene_type):

        self.evaluate(population, eval_func)
        parent1, parent2 = self.__select_parents(population)
        child1, child2 = self.__crossover(
            parent1['gene'], parent2['gene'], lethal_gene,
            death_func)
        if random.randint(0, 1) == 0:
            child1 = self.__mutate(
                child1['gene'], lethal_gene, death_func, g_min,
                g_max, gene_type)
        else:
            child2 = self.__mutate(
                child2['gene'], lethal_gene, death_func, g_min,
                g_max, gene_type)

        family = [parent1, parent2, child1, child2]
        self.evaluate(family, eval_func)
        score = self.get_score(family)

        ic, sc = min(score[2:]), max(score[2:])
        ip, sp = min(score[:-2]), max(score[:-2])
        if ic >= sp:
            population.extend([child1, child2])
            population.remove(family[score.index(ip)])
        elif sc < ip:
            population.remove(family[score.index(ip)])
        elif sc < sp:
            population.append(family[score.index(sc)])
            population.remove(family[score.index(ip)])
        elif sc >= sp:
            population.append(family[score.index(sc)])
            population.remove(parent1)
            population.remove(parent2)
            population.append(
                {'score': None, 'gene': self.__get_population(
                    1, n_gene, g_min, g_max, gene_type)})

        if len(population) == 1:
            population.append(
                {'score': None, 'gene': self.__get_population(
                    1, n_gene, g_min, g_max, gene_type)})

        self.evaluate(population, eval_func)

        return population

    def PfGA_main(
            self, n_generation, n_gene, eval_func, lethal_gene,
            death_func, g_min, g_max, gene_type):

        population = [
            {'score': None, 'gene': p}
            for p in self.__get_population(
                2, n_gene, g_min, g_max, gene_type)
        ]

        for g in range(n_generation):
            print('Generation: ' + str(g + 1))

            population = self.generation_change(
                population, eval_func, n_gene, lethal_gene,
                death_func, g_min, g_max, gene_type)

            score = self.get_score(population)

            print('min: ', min(score))  # 最も低い評価を表示
            print('max: ', max(score))  # 最も高い評価を表示
            print('ave: ', statistics.mean(score))  # 評価の平均値を表示
            print()

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)
        print('best')
        print(population[0])  # 最も高い評価の個体を表示


class PfGA_binary(PfGA, PfGA_method):

    def __init__(self,
                 n_gene,
                 n_generation,
                 eval_func=sum,
                 lethal_gene=False,
                 death_func=None):

        super().__init__(
            n_gene, n_generation, eval_func, lethal_gene,
            death_func)

    def main(self):

        super().PfGA_main(
            self.N_GENERATION, self.N_GENE,
            self.EVALUATE_FUNCTION, self.LETHAL_GENE,
            self.DEATH_FUNCTION, _, _, 'binary')


class PfGA_real(PfGA, PfGA_method):

    def __init__(self,
                 n_gene,
                 n_generation,
                 eval_func,
                 g_min,
                 g_max,
                 lethal_gene=False,
                 death_func=None):

        super().__init__(
            n_gene, n_generation, eval_func, lethal_gene,
            death_func)
        self.__GENE_MIN = g_min  # 遺伝子の最小値
        self.__GENE_MAX = g_max  # 遺伝子の最大値

    def main(self):

        super().PfGA_main(
            self.N_GENERATION, self.N_GENE,
            self.EVALUATE_FUNCTION, self.LETHAL_GENE,
            self.DEATH_FUNCTION, self.__GENE_MIN, self.__GENE_MAX,
            'real')


if __name__ == '__main__':
    GA_binary(10, 20, 25).main()
    # GA_island_binary(10, 20, 25, 5, 0.2, 4).main()
    # PfGA_binary(10, 25).main()

# http://testpy.hatenablog.com/entry/2017/01/09/213846
# http://qiita.com/Azunyan1111/items/975c67129d99de33dc21
# http://ichitcltk.hustle.ne.jp/gudon2/index.php?pageType=file&id=python_class_inheritance.md
# http://www.sist.ac.jp/~kanakubo/research/evolutionary_computing/parameter_free_ga.html
# http://qiita.com/simanezumi1989/items/4f821de2b77850fcf508