import numpy as np
import pandas as pd
import math
import random

def startTrain(PredictClass, TrueClass):
    n = PredictClass.shape[1]
    result = PEP(n, PredictClass, TrueClass)
    return result


def objFunc(solution, PredictClass, TrueClass):
    subPredInd = np.where(solution == 1)[0]
    subLength = len(subPredInd)
    subMat = PredictClass[:, subPredInd]
    res = np.zeros((PredictClass.shape[0], subLength))

    for i in range(0, subLength):
        res[:, i] = np.where(subMat[:, i] == TrueClass, 1, -1)
    result = np.sum(res, axis=1)

    fValue = np.sum(result < 0) + np.sum(result == 0) / 2
    return fValue

def eval(solution, PredictClass, TrueClass):
    subPredInd = np.where(solution == 1)[0]
    subLength = len(subPredInd)
    subMat = PredictClass[:, subPredInd]
    res = np.zeros((PredictClass.shape[0], subLength))

    for i in range(0, subLength):
        res[:, i] = np.where(subMat[:, i] == TrueClass, 1, -1)
    result = np.sum(res, axis=1)

    evalValue = np.sum(result < 0) + np.sum(result == 0) / 2
    return evalValue


def PEP(n, PredictClass, TrueClass):
    population = pd.DataFrame(np.random.randint(0, 2, n)) # integers in [0,2)
    popSize = 1
    fitness = pd.DataFrame(np.zeros(2))
    fitness.iloc[1,0] = sum(population.iloc[:,0].values)
    if fitness.iloc[1, 0] == 0:
        fitness.iloc[0, 0] = float('inf')
    else:
        fitness.iloc[0, 0] = objFunc(population, PredictClass, TrueClass)

    T = round((n ** 2) * math.log(n))
    for it in range(0, T):
        offspring = abs(population.iloc[:,np.random.randint(0, popSize)].values - np.random.choice([1,0],size=n,p=[1/n,1-1/n]))
        offspringFit = np.zeros(2)
        offspringFit[1] = sum(offspring)
        if offspringFit[1] == 0:
            offspringFit[0] = float('inf')
        else:
            offspringFit[0] = objFunc(offspring, PredictClass, TrueClass)
        domain_count = 0
        deleteIndex = []
        for i in range(0, popSize):
            if fitness.iloc[0, i] < offspringFit[0] and fitness.iloc[1, i] <= offspringFit[1]:
                domain_count = domain_count + 1
            if fitness.iloc[0, i] <= offspringFit[0] and fitness.iloc[1, i] < offspringFit[1]:
                domain_count = domain_count + 1

        if domain_count > 0:
            continue
        else:
            for j in range(0, popSize):
                if fitness.iloc[0, j] >= offspringFit[0] and fitness.iloc[1, j] >= offspringFit[1]:
                    deleteIndex.append(j)

        population.drop(population.columns[deleteIndex], axis=1, inplace=True)
        fitness.drop(fitness.columns[deleteIndex], axis=1, inplace=True)
        popSize = population.shape[1]
        colname = list(range(popSize))
        population.colnames = colname
        fitness.colnames = colname
        population.loc[:, popSize] = offspring
        fitness.loc[:, popSize] = offspringFit
        popSize = population.shape[1]

        bestSolution = np.zeros(n)
        bestFitness = np.array([float('inf'), float('inf')])
        isChanged = np.zeros(n)
        for j in range(0, n):
            hammingSolutions = np.tile(offspring, (n - j, 1))
            hammingFValues = np.zeros(n - j)
            temp = np.where(isChanged == 0)[0]
            for p in range(0, n-j):
                hammingSolutions[p, temp[p]] = 1 - hammingSolutions[p, temp[p]]
                hammingFValues[p] = objFunc(hammingSolutions[p], PredictClass, TrueClass)
                if 0 == sum(hammingSolutions[p]):
                    hammingFValues[p] = float('inf')

            bestPos = np.where(hammingFValues == np.min(hammingFValues))[0][0]
            bestF = hammingFValues[bestPos]
            bestS = sum(hammingSolutions[bestPos])


            if bestF < bestFitness[0]:
                bestSolution = hammingSolutions[bestPos]
                bestFitness = [bestF, bestS]
            else:
                if bestF <= bestFitness[0] and bestS <= bestFitness[1]:
                    bestSolution = hammingSolutions[bestPos]
                    bestFitness = [bestF, bestS]
            isChanged[np.where(offspring != hammingSolutions[bestPos])[0]] = 1
            offspring = hammingSolutions[bestPos]
        dc = 0
        deleteInd = []
        for i in range(0, popSize):
            if fitness.iloc[0, i] < bestFitness[0] and fitness.iloc[1, i] <= bestFitness[1]:
                dc = dc + 1
            if fitness.iloc[0, i] <= bestFitness[0] and fitness.iloc[1, i] < bestFitness[1]:
                dc = dc + 1

        if dc > 0:
            continue
        else:
            for k in range(0, popSize):
                if fitness.iloc[0, k] >= bestFitness[0] and fitness.iloc[1, k] >= bestFitness[1]:
                    deleteInd.append(k)

        population.drop(population.columns[deleteInd], axis=1, inplace=True)
        fitness.drop(fitness.columns[deleteInd], axis=1, inplace=True)
        popSize = population.shape[1]
        colname = list(range(popSize))
        population.colnames = colname
        fitness.colnames = colname
        population.loc[:, popSize] = bestSolution
        fitness.loc[:, popSize] = bestFitness
        popSize = population.shape[1]

    evalValue = np.zeros(popSize)
    for i in range(0, popSize):
        evalValue[i] = eval(population.iloc[:,i], PredictClass, TrueClass)

    index = np.where(evalValue == np.min(evalValue))[0][0]
    selectedEnsemble = population.iloc[:,index]
    return selectedEnsemble