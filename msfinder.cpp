#include<iostream>
#include<cstdlib>
#include<random>
#include<algorithm>
using namespace std;

class MagicSquares{
    private:
        int m_N;
        int m_arrN;
        int m_size;
        int** m_squares;
        int* m_fitness;
    public:
        MagicSquares(int N, int size);
        void calFit();
        void printFirst();
        int getFit(int pos);
        void Sort();
        void breed();
        void mutate(int i);
        void copy(int i, int j);
        void swap(int i, int j);
        void crossOver(int n, int i, int j);
        void print(int i);
};

int main(){

    //reading N
    int N=5;
    int gens=1000;

    //Initialize a population of magic squares with randomly generated values
    MagicSquares squares(N,100);

    long int i=0;

    //for(int i = 0; i < gens; ++i){
    while(1){
        //Evaluate the fitness of each square in the population, 
        //using a function that measures how close the square is to being a true magic square
        squares.calFit();

        /*for(int i = 0; i<20; ++i){
            cout << i << ". " << squares.getFit(i) << endl;
            squares.print(i);
            cout << endl;
        }*/

        //Select the most fit individuals from the current population to serve as parents for the next generation
        //Breed the selected parents to create a new population of magic squares, 
        //using crossover and mutation to introduce genetic variation.
        squares.breed();

        //cout << "Iteration: " << i << " Fit: " << squares.getFit(0) << " _ ";
        //squares.printFirst();
        //cout << endl;

        if(i%1000==0){
            cout << i << ":  fitness: " << squares.getFit(0) << endl;
        }
        ++i;
    }
    return 0;
}

MagicSquares::MagicSquares(int N, int size){//creating random squares
    m_N = N;
    m_size = size;
    m_arrN = m_N*m_N;
    m_squares = new int*[m_size];
    m_fitness = new int[m_size];
    random_device rd;
    mt19937 gen(rd());
    for (int i = 0; i < m_size; i++) {
        m_squares[i] = new int[m_arrN];
        for (int j = 0; j < m_arrN; j++) {
            m_squares[i][j] = j+1;
        }
        shuffle(m_squares[i], m_squares[i]+m_arrN, gen);
    }

    /*for(int i=0; i<size; ++i){
        for(int j=0; j<m_arrN; ++j){
            cout << m_squares[i][j] << " ";
        }
        cout << endl;
    }*/
}

void MagicSquares::calFit(){//calculate the fitness of all squares
    int error;
    int first_row;
    int val;
    for(int i = 0; i < m_size; ++i){
        first_row = 0;
        //first row as set point
        for(int j = 0; j < m_N; ++j){
            first_row += m_squares[i][j];
        }
        error = 0;
        //error in rows
        for(int j = 1; j < m_N; ++j){
            val = 0;
            for(int k = 0; k < m_N; ++k){
                val += m_squares[i][j*m_N+k];
            }
            error += abs(first_row-val);
        }
        //error in columns
        for(int j = 0; j < m_N; ++j){
            val = 0;
            for(int k = 0; k < m_N; ++k){
                val += m_squares[i][k*m_N+j];
            }
            error += abs(first_row-val);
        }
        m_fitness[i]=error;
        if(error==0){
            this->print(i); //if error=0 the square is magic, print it!
        }
    }
}

void MagicSquares::printFirst(){
    for(int j=0; j<m_arrN; ++j){
        cout << m_squares[0][j] << " ";
    }
}

int MagicSquares::getFit(int pos){
    return m_fitness[pos];
}

void MagicSquares::Sort(){
    for(int i = 0; i < m_size; ++i){
        for(int j = 0; j< m_size-i-1; ++j){
            if(m_fitness[j]>m_fitness[j+1]){
                this->swap(j,j+1);
            }
        }
    }
}

void MagicSquares::mutate(int i){//mutate i by swapping two random numbers
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, m_arrN-1);
    int index1 = dis(gen);
    int index2 = dis(gen);
    int buffer = m_squares[i][index1];
    m_squares[i][index1] = m_squares[i][index2];
    m_squares[i][index2] = buffer;
}

void MagicSquares::copy(int i, int j){ //from i to j (j gets ereased)
    for(int k = 0; k < m_arrN; ++k){
        m_squares[j][k] = m_squares[i][k];
    }
}

void MagicSquares::swap(int i, int j){
    int buffer;
    for(int k = 0; k < m_arrN; ++k){
        buffer = m_squares[i][k];
        m_squares[i][k] = m_squares[j][k];
        m_squares[j][k] = buffer;
    }
    int buf = m_fitness[i];
    m_fitness[i] = m_fitness[j];
    m_fitness[j] = buf;
}

void MagicSquares::breed(){
    /*for(int i=0; i<m_size; ++i){
        cout << m_fitness[i] << " ";
    }
    cout << endl;*/
    this->Sort();
    /*for(int i=0; i<m_size; ++i){
        cout << m_fitness[i] << " ";
    }
    cout << endl;*/
    //fill the later (worse) squares with either mutation of "good ones"
    //or crossover of two "good" ones
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);
    uniform_int_distribution<> dis1(0, m_size/4-1);
    for(int i = m_size/4; i<m_size; ++i){
        //cout << "breed1, " << i << ": "; 
        if(dis(gen)==0){ //mutation of a "good" square (first copy then mutate)
            //cout << "mutate" <<endl;
            this->copy(dis1(gen),i);
            this->mutate(i);
        }
        else{
        //cout << "crossover" << endl;
            this->crossOver(i,dis1(gen),dis1(gen)); //crossover of two "good" squares
        }
    }

    //mutate the best 25% (saved first in the array)
    for(int i = 0; i<m_size/4; ++i){
        //cout << "breed 2, " << i << endl;
        this->mutate(i);
    }
}

void MagicSquares::crossOver(int n, int i, int j){//save crossover of i and j in n
    int invi[m_arrN];
    int invj[m_arrN];
    int invn[m_arrN];

    for(int k = 0; k < m_arrN; ++k){
        invi[k]=0;
        invj[k]=0;
    }

    //calculating inversion sequences of i and j
    for(int k1=1; k1<=m_arrN; ++k1){
        for(int k2=0; m_squares[i][k2]!=k1; ++k2){
            if(m_squares[i][k2]>k1){
                invi[k1]++;
            }
        }
        for(int k2=0; m_squares[j][k2]!=k1; ++k2){
            if(m_squares[j][k2]>k1){
                invj[k1]++;
            }
        }
    }

    //generating random pivot index
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, m_arrN-1);
    int pivot = dis(gen);

    for(int k = 0; k < pivot; ++k){
        invn[k]=invi[k];
    }
    for(int k = pivot; k < m_arrN; ++k){
        invn[k]=invj[k];
    }

    //generating n with its inversion sequence
    for(int k = 0; k < m_arrN; ++k){
        m_squares[n][k]=0;
    }
    for(int k1=1; k1<=m_arrN; ++k1){
        for(int k2=0; k2<m_arrN; ++k2){
            if(m_squares[n][k2]==0 && invn[k1-1]==0){
                m_squares[n][k2]=k1;
                break;
            }
            else{
                if(m_squares[n][k2]==0){
                    invn[k1-1]--;
                }
                else{
                    if(invn[k1-1]!=0){
                        invn[k1-1]--;
                    }
                }
            }
        }
    }

}

void MagicSquares::print(int i){
    cout << "Magic Square Alarm!" << endl;
    for(int k = 0; k < m_N; ++k){
        for(int j = 0; j < m_N; ++j){
            cout << m_squares[i][k*m_N+j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    exit(0);
}