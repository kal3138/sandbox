# -*- coding: utf-8 -*-

import statistics

import de


class eDE_params(DE_params):

    def __init__(self,
                 n_gene,
                 n_po,
                 n_generation,
                 g_min,
                 g_max,
                 eval_func,
                 rest_func):

        super().__init__(eval_func,
                         n_gene,
                         n_po,
                         n_generation,
                         g_min,
                         g_max,
                         eval_func)
        self.RESTRICT_FUNCTION = rest_func  # 制約関数


class eDE_evaluate(DE_evaluate):

    def get_violation(self, population):

        violation_list = [None] * len(population)
        for i, p in enumerate(population):
            violation_list[i] = p['violation']

        return violation_list

    def calc_violation(self, x, rest_func):

        violation = 0
        for rf in rest_func:
            violation += rf(x)

        return violation

    def evaluate(self, population, eval_val, rest_func):

        for p in population:
            if p['score'] is None:
                p['score'] = super().calc_score(
                    p['gene'], eval_val)
            if p['violation'] is None:
                p['violation'] = self.calc_violation(
                    p['gene'], rest_func)


class eDE_method(DE_method, eDE_evaluate):

    def update_epsilon(self, generation, n_generation):

        pass

        return eps

    def eDE_main(self, n_generation, n_pop, n_gene, gene,
                 eval_func, rest_func, f, cr, crs, base, p, f2):
        population = [
            {'score': None, 'violation': None, 'gene': p}
            for p in super().get_population(n_pop, n_gene, gene)
        ]
        super().evaluate(population, eval_func, rest_func)

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
                violation_child = super().calc_violation(
                    child, rest_func)
                if (violation_child > eps) ^ \
                        (parent['violation'] > eps):
                    population.append(
                        parent if parent['violation'] <= eps else
                        {
                            'score': score_child,
                            'violation': violation_child,
                            'gene': child})
                else:
                    if score_child > parent[0]:
                        population.append((score_child, child))
                    else:
                        population.append(parent)

            score = [
                s for s, v in zip(
                    super().get_score(population),
                    super().get_violation(population)
                ) if v <= eps]

            print('min: ', min(score))  # 最も低い評価を表示
            print('max: ', max(score))  # 最も高い評価を表示
            print('ave: ', statistics.mean(score))  # 評価の平均値を表示

            eps = update_epsilon(g, n_generation)

        population.sort(
            key=lambda x: (x['score'] is not None, x),
            reverse=True)

        return population


class eDE(eDE_params, eDE_method):

    def __init__(self, n_gene, n_po, n_generation, g_min, g_max,
                 eval_func, rest_func, f, cr):

        super().__init__(n_gene, n_po, n_generation, g_min, g_max,
                         eval_func, rest_func, f, cr)

    def main(self):

        gene = (self.GENE_MIN, self.GENE_MAX)
        population = super().eDE_main(
            self.N_GENERATION, self.N_POP, self.N_GENE, gene,
            self.EVALUATE_FUNCTION, self.RESTRICT_FUNCTION,
            self.SCALING_PARAMETER, self.CROSSOVER_RATE,
            self.CROSSOVER_TYPE, self.BASE, self.TOP_PERCENT,
            self.SCALING_PARAMETER2)

        print('best')
        print(population[0])  # 最も高い評価の個体を表示

# http://www.ints.info.hiroshima-cu.ac.jp/~kushida/ML/ML11.pdf
# https://qiita.com/LeftLetter/items/124461ff2f0841d9132b