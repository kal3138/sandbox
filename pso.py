# -*- coding: utf-8 -*-

import random
import copy


class PSO_params(object):

    value_LDIWM = {
        'cr': 2.0,
        'cs': 2.0
    }
    value_CM = {
        'w': lambda g: 0.729,
        'cr': 1.4955,
        'cs': 1.4955
    }
    value_RIWM = {
        'w': lambda g: random.uniform(0.5, 1),
        'cr': 1.4955,
        'cs': 1.4955
    }

    def __init__(
            self, n_dim, n_po, n_generation, p_min, p_max, v_max,
            eval_func):

        self.N_DIM = n_dim  # 遺伝子長
        self.N_POP = n_po  # 1世代の個体数
        self.N_GENERATION = n_generation  # 最大計算世代数
        self.POSITION_MIN = p_min  # 位置の最小値
        self.POSITION_MAX = p_max  # 位置の最大値
        self.VELOCITY_MAX = v_max  # 速度の最大値
        self.EVALUATE_FUNCTION = eval_func  # 評価関数

        # set_paramで変更可能
        # PSOのパラメータ設定
        self.INERTIA_WEIGHT = lambda g: 0.9 - 0.5 / self.N_GENERATION * g
        self.RECOGNISION_PARAMETER = value_LDIWM['cr']
        self.SOCIETY_PARAMETER = value_LDIWM['cs']
        self.SYNC_RANDOM = False
        self.PSO_BEST = 'gbest'
        self.NEIGHBORHOOD_SIZE = 1

        try:
            super().__init__()
        except AttributeError:
            pass

    def set_param(
            self,
            param_mode=None,
            w=None,
            cr=None,
            cs=None,
            sr=None,
            pso_best=None,
            neis=1):

        if param_mode is not None:
            if param_mode == 'LDIWM':
                pass

            elif param_mode == 'CM':
                self.INERTIA_WEIGHT = value_CM['w']
                self.RECOGNISION_PARAMETER = value_CM['cr']
                self.SOCIETY_PARAMETER = value_CM['cs']

            elif param_mode == 'RIWM':
                self.INERTIA_WEIGHT = value_RIWM['w']
                self.RECOGNISION_PARAMETER = value_RIWM['cr']
                self.SOCIETY_PARAMETER = value_RIWM['cs']

            else:
                print('Unknown Mode')
                return

        else:
            if w is not None:
                self.INERTIA_WEIGHT = w
            if cr is not None:
                self.RECOGNISION_PARAMETER = cr
            if cs is not None:
                self.SOCIETY_PARAMETER = cs

        if sr is not None:
            self.SYNC_RANDOM = sr

        if pso_best is not None:
            if pso_best == 'GBEST':
                pass

            elif pso_best == 'LBEST':
                self.NEIGHBORHOOD_SIZE = neis

            elif pso_best == 'INSM':
                self.NEIGHBORHOOD_SIZE = 3

            else:
                print('Unknown Mode')
                return

        try:
            super().set_param(
                param_mode, w, cr, cs, sr, pso_best, neis)
        except AttributeError:
            pass


class PSO_evaluate(object):

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

    def evaluate(self, particle, eval_val):

        if particle['score'] is None:
            particle['score'] = self.calc_score(
                particle['position'], eval_val)
        if particle['score_best'] is None:
            particle['score_best'] = particle['score']
            particle['position_best'] = particle['position']

        try:
            super().evaluate(particle, eval_val)
        except AttributeError:
            pass

    def evaluate_best(self, particle, p_best):

        if particle['score'] > particle['score_best']:
            particle['score_best'] = particle['score']
            particle['position_best'] = particle['position']
            if particle['score_best'] > p_best['score_best']:
                p_best = particle

        try:
            super().evaluate_best(particle, p_best)
        except AttributeError:
            pass


class PSO_method(PSO_evaluate):

    def __get_population(self, n_pop, n_dim, pos, v_max):

        ru = random.uniform
        p_min, p_max = pos
        v_min, v_max = vel
        for i in range(n_pop):
            position, velocity = \
                [ru(p_min, p_max) for _ in range(n_dim)], \
                [ru(v_min, v_max) for _ in range(n_dim)]

            yield position, velocity

    def __update_velocity(
            self, n_dim, particle, g_pos_best, w, cr, cs, v_max,
            sync_rand):

        rr = random.random
        if sync_rand:
            r1 = rr()
            r2 = rr()
        else:
            r1 = [rr() for _ in range(n_dim)]
            r2 = [rr() for _ in range(n_dim)]

        v_update = w * particle['velocity'] \
            + cr * r1 * (
                particle['position_best'] -
                particle['position']) \
            + cs * r2 * (g_pos_best - particle['position'])
        if v_update < -v_max:
            v_update = -v_max
        elif v_update > v_max:
            v_update = v_max

        particle['velocity'] = v_update

    def __update_position(self, particle):

        particle['position'] += particle['velocity']

    def __init_gbest(self, population):

        return max(
            population, key=lambda x: (x['score'] is not None, x))

    def __init_lbest(self, population, neighborhood_size):

        n = len(population)
        j = neighborhood_size // 2
        l_best = [None] * n
        for i in range(n):
            k = i if i + j < n else i - n
            neighborhood = population[k - j:k + j + 1]
            l_best[i] = max(
                neighborhood,
                key=lambda x: (x['score'] is not None, x))

        return l_best

    def __select_pbest(self, population, neighborhood_size=1):

        if neighborhood_size == 1:
            p_best = self.__init_gbest(population)
        elif neighborhood_size > 1:
            p_best = self.__init_lbest(population, neighborhood_size)

        return p_best

    def __update_evaluate(
            self, n_dim, particle, p_best, eval_func, w, cr, cs,
            v_max, sync_rand):

        self.__update_velocity(
            n_dim, particle, p_best['position_best'], w(g), cr,
            cs, v_max, sync_rand)
        self.__update_position(p)
        super().evaluate(particle, eval_func)
        super().evaluate_best(particle, p_best)

    def PSO_main(
            self, n_generation, n_pop, n_dim, pos, v_max,
            eval_func, w, cr, cs, sync_rand, neighborhood_size):
        population = [{
            'position': p,
            'velocity': v,
            'score': None,
            'score_best': None,
            'position_best': [None] * n_dim}
            for p, v in super().__get_population(
                n_pop, n_dim, pos, v_max)
        ]
        for p in range(population):
            super().evaluate(p, eval_func)

        p_best_list = self.__select_pbest(
            population, neighborhood_size)

        if len(p_best_list) == 1:
            p_best = p_best_list
            for g in range(n_generation):
                print('Generation: ' + str(g + 1))
                for p in population:
                    self.__update_evaluate(
                        n_dim, p, p_best, eval_func, w, cr, cs,
                        v_max, sync_rand)
        else:
            for g in range(n_generation):
                print('Generation: ' + str(g + 1))
                for i, p in enumerate(population):
                    p_best = p_best_list[i]
                    self.__update_evaluate(
                        n_dim, p, p_best, eval_func, w, cr, cs,
                        v_max, sync_rand)
                    p_best_list[i] = p_best

            print('best: ', max(
                p_best_list,
                key=lambda x: (x['score_best'] is not None, x)))    # 最も高い評価を表示

        population.sort(
            key=lambda x: (x['score_best'] is not None, x),
            reverse=True)

        return population


class PSO(PSO_params, PSO_method):
    def __init__(self,
                 n_dim,
                 n_po,
                 n_generation,
                 p_min,
                 p_max,
                 v_max,
                 eval_func):
        super().__init__(
            n_dim, n_po, n_generation, p_min, p_max, v_max,
            eval_func)

    def main(self):
        pos = (self.POSITION_MIN, self.POSITION_MAX)
        population = super().PSO_main(
            self.N_GENERATION, self.N_POP, self.N_DIM, pos,
            self.VELOCITY_MAX, self.EVALUATE_FUNCTION,
            self.INERTIA_WEIGHT, self.RECOGNISION_PARAMETER,
            self.SOCIETY_PARAMETER, self.SYNC_RANDOM,
            self.NEIGHBORHOOD_SIZE)

        print('best')
        print(population[0])  # 最も高い評価の個体を表示

# http://www.ints.info.hiroshima-cu.ac.jp/~kushida/ML/ML12.pdf