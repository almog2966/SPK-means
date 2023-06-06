import sys
import pandas as pd
import numpy as np
import myspkmeans
import kmeans_c

E = 0
MAX_ITER = 300

#kmeans func from HW2
def kmeans_pp(k,dataPoints):
    n_rows,n_col = dataPoints.shape
    if(k>= n_rows):
        print("An Error Has Occurred")
        exit(1)
    np.random.seed(0)
    index_Arary = []
    centroids = np.zeros([k, dataPoints.shape[1]])
    centerIndex = np.random.choice(n_rows)
    index_Arary.append(centerIndex)
    centroids[0] = dataPoints[centerIndex]
    i = 1
    while i<k:
        D_x = np.sqrt(np.sum(np.power(
            np.swapaxes(dataPoints * np.ones([i, *dataPoints.shape]), 0, 1)
            -
            centroids[:i]
            , 2), axis=2))
        c = centroids[:i]
        D_x = np.min(D_x,axis = 1)
        prob_Array = D_x/np.sum(D_x)
        centerIndex = np.random.choice(a = dataPoints.shape[0], p = prob_Array)
        index_Arary.append(centerIndex)
        centroids[i] = dataPoints[centerIndex]
        i+=1
    return centroids, index_Arary

#reading inputs from args array then reading file
#return file as numpy list, k if got one and goal
def read_inputs():
    if(len(sys.argv) ==4):
        k,goal,file_name = int(sys.argv[1]) , sys.argv[2], sys.argv[3]
    elif(len(sys.argv) == 3):
        goal,file_name = sys.argv[1] , sys.argv[2]
        k = -1
    else:
        return -1,-1,-1
    try:
        file_data = pd.read_csv(file_name , sep = ',', header=None)
    except:
        print("An Error Has Occurred")
        exit(1)
    return k , goal, file_data.values.tolist()


#print array each value as .4f
def printVector(vector):
    for i in range(len(vector)):
        if(i != len(vector)-1):
            print(format(vector[i],'.4f'),",",end='',sep='')
        else:
            print(format(vector[i],'.4f'))

#print int array 
def printIndex(vector):
    for i in range(len(vector)):
        if(i != len(vector)-1):
            print(vector[i],",",end='',sep='')
        else:
            print(vector[i])

#print first row as eigenvalues other row is eigen vector as matrix
def printJacobi(matrix):
    printVector(matrix[len(matrix)-1])
    for i in range(len(matrix)-1):
        for j in range(len(matrix[i])):
            if j == len(matrix[i])-1:
                print(format(matrix[i][j],'.4f'))
            else:
                print(format(matrix[i][j],'.4f'),",",end='',sep='')

#print 2D array each value as .4f
def printMatrix(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if j == len(matrix[i])-1:
                print(format(matrix[i][j],'.4f'))
            else:
                print(format(matrix[i][j],'.4f'),",",end='',sep='')
       



def main():
    k, goal, matrix = read_inputs()
    if(matrix == -1):
        print("An Error Has Occurred")
        return 0
    if goal == "spk":
        dataPoints = myspkmeans.spk(k,matrix)

        if(k==-1):
            k = len(dataPoints[0])
        n = len(dataPoints)
        centroids,index_Arary = kmeans_pp(k,np.array(dataPoints))
        centroids_ret = kmeans_c.fit(k, n, k, MAX_ITER, E, dataPoints, centroids.tolist())
        printIndex(index_Arary)
        printMatrix(centroids_ret)
    elif goal == "wam":
        retMatrix = myspkmeans.wam(matrix)
        printMatrix(retMatrix)
    elif goal == "ddg":
        retMatrix = myspkmeans.ddg(matrix)
        printMatrix(retMatrix)
    elif goal == "gl":
        retMatrix = myspkmeans.gl(matrix)
        printMatrix(retMatrix)
    elif goal == "jacobi":
        retMatrix = myspkmeans.jacobi(matrix)
        printJacobi(retMatrix)
    else:
        print("An Error Has Occurred")


if __name__ == '__main__':
    main()