# -*- coding: utf-8 -*-

import random
import math
import statistics
import copy


class DE_params(object):

    default_value = {
        'f': 0.5,
        'cr': 0.5,
        'crs': 'bin',
        'base': 'rand1',
        'p': None,
        'f2': None
    }

    def __init__(self,
                 n_gene,
                 n_po,
                 n_generation,
                 g_min,
                 g_max,
                 eval_func):

        self.N_GENE = n_gene  # 遺伝子長
        self.N_POP = n_po  # 1世代の個体数
        self.N_GENERATION = n_generation  # 最大計算世代数
        self.GENE_MIN = g_min  # 遺伝子の最小値
        self.GENE_MAX = g_max  # 遺伝子の最大値
        self.EVALUATE_FUNCTION = eval_func  # 評価関数

        # set_paramで変更可能
        # スケーリングパラメータ
        self.SCALING_PARAMETER = DE_params.default_value['f']
        # 交叉率
        self.CROSSOVER_RATE = DE_params.default_value['cr']
        # 交叉方法
        self.CROSSOVER_TYPE = DE_params.default_value['crs']
        # ベースベクトル選択方法
        self.BASE = DE_params.default_value['base']
        # 上位層のパーセンテージ(current-to-p_best)
        self.TOP_PERCENT = DE_params.default_value['ｐ']
        # スケーリングパラメータ
        self.SCALING_PARAMETER2 = DE_params.default_value['f2']

        try:
            super().__init__()
        except AttributeError:
            pass

    def set_param(
            self,
            f=None,
            cr=None,
            crs=None,
            base=None,
            p=None,
            f2=None,
    ):

        if f is not None:
            self.SCALING_PARAMETER = f
        if mu is not None:
            self.CROSSOVER_RATE = cr
        if el is not None:
            self.CROSSOVER_TYPE = crs
        if nt is not None:
            self.BASE = base
        if pe is not None:
            self.TOP_PERCENT = p
        if lg is not None:
            self.SCALING_PARAMETER2 = f2

        try:
            super().set_param(f, cr, crs, base, p, f2)
        except AttributeError:
            pass


class DE_evaluate(object):

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

    def evaluate(self, population, eval_val):

        for p in population:
            if p['score'] is None:
                p['score'] = self.calc_score(p['gene'], eval_func)

        try:
            super().evaluate(population)
        except AttributeError:
            pass

        return fitness


class DE_method(DE_evaluate):

    def __get_population(self, n_pop, n_gene, gene):

        g_min, g_max = gene
        for i in range(n_pop):
            population = [
                random.uniform(g_min, g_max)
                for _ in range(n_gene)]

            yield population

    def __select_method(self, population, base, trg, p):

        if base == 'rand1':
            selected = random.sample(population, 3)
        elif base == 'rand2':
            selected = random.sample(population, 5)
        elif base == 'best':
            population.sort(
                key=lambda x: (x['score'] is not None, x),
                reverse=True)
            selected = [None] * 3
            selected[0] = population[0]
            selected[1:] = random.sample(population[1:], 2)
        elif base == 'current-to':
            if trg is None:
                print('target is None')
                return
            selected = [None] * 3
            selected[0] = trg
            selected[1:] = random.sample(population, 2)
        elif any(base == b for b in (
                'current-to-best', 'current-to-p_best')):
            if trg is None:
                print('target is None')
                return
            population.sort(
                key=lambda x: (x['score'] is not None, x),
                reverse=True)
            selected = [None] * 5
            selected[0:3:2] = trg
            if base == 'current-to-best':
                selected['gene'] = population[0]
                selected[3:] = random.sample(population[1:], 2)
            else:
                s = random.sample(
                    population[
                        :math.ceil(len(population) * p / 100)], 3)
                s.sort(
                    key=lambda x: (x['score'] is not None, x),
                    reverse=True)
                selected['gene'] = s[0]
                selected[3:] = s[1:]
        else:
            print('Unknown base')
            return

        return selected

    def __mutate(self, selected, f, base, f2):

        if any(base == b for b in (
                'rand1', 'best', 'current-to')):
            parent1, parent2, parent3 = selected
            mutant = parent1['gene'] + \
                f * (parent2['gene'] - parent3['gene'])
        elif base == 'rand2':
            parent1, parent2, parent3, parent4, parent5 = selected
            mutant = parent1['gene'] + \
                f * (parent2['gene'] - parent3['gene']) + \
                f * (parent4['gene'] - parent5['gene'])
        elif any(base == b for b in (
                'current-to-best', 'current-to-p_best')):
            if f2 is None:
                print('Scaling Factor is None')
                return
            parent1, parent2, parent3, parent4, parent5 = selected
            mutant = parent1['gene'] + \
                f2 * (parent2['gene'] - parent3['gene']) + \
                f * (parent4['gene'] - parent5['gene'])
        else:
            print('Unknown base')
            return

        return mutant

    def __crossover(self, parent, mutant, cr, crs):

        l_parent = len(parent)
        rd = random.random
        j = random.randint(0, l_parent - 1)
        if crs == 'bin':
            child = [
                gm if rd() < cr or i == j else g1
                for i, (gp, gm) in enumerate(zip(parent, mutant))
            ]
        elif crs == 'exp':
            child = copy.deepcopy(parent)
            for k in range(j + 1, l_parent):
                if rd() >= cr:
                    child[j:k - 1] = mutant[j:k - 1]
                    break
            else:
                child[j:] = mutant[j:]
        else:
            print('Unknown mode')
            return

        return child

    def DE_main(
            self, n_generation, n_pop, n_gene, gene, eval_func, f,
            cr, crs, base, p, f2):

        population = [
            {'score': None, 'gene': p}
            for p in super().get_population(n_pop, n_gene, gene)
        ]
        super().evaluate(population, eval_func)

        for g in range(n_generation):
            print('Generation: ' + str(g + 1))
            for i in range(n_pop):
                parent = population[i]
                population.pop(i)

                selected = self.__select_method(
                    population, base, parent, p)
                mutant = self.__mutate(selected, f, base, f2)
                child = self.__crossover(
                    parent['gene'], mutant, cr, crs)
                score_child = super().calc_score(child, eval_func)
                if score_child > parent[0]:
                    population.append((score_child, child))
                else:
                    population.append(parent)

            score = super().get_score(population)

            print('min: ', min(score))  # 最も低い評価を表示
            print('max: ', max(score))  # 最も高い評価を表示
            print('ave: ', statistics.mean(score))  # 評価の平均値を表示

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)

        return population


class DE(DE_params, DE_method):

    def __init__(
            self, n_gene, n_po, n_generation, g_min, g_max,
            eval_func, f, cr):

        super().__init__(
            n_gene, n_po, n_generation, g_min, g_max, eval_func,
            f, cr)

    def main(self):

        gene = (self.GENE_MIN, self.GENE_MAX)
        population = super().DE_main(
            self.N_GENERATION, self.N_POP, self.N_GENE, gene,
            self.EVALUATE_FUNCTION, self.SCALING_PARAMETER,
            self.CROSSOVER_RATE, self.CROSSOVER_TYPE, self.BASE,
            self.TOP_PERCENT, self.SCALING_PARAMETER2)

        print('best')
        print(population[0])  # 最も高い評価の個体を表示

# http://www.ints.info.hiroshima-cu.ac.jp/~kushida/ML/ML11.pdf
# https://qiita.com/LeftLetter/items/124461ff2f0841d9132b