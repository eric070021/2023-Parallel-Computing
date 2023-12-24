#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <mpi.h>

using namespace std;

#define MIN(a,b)   ((a)<(b)?(a):(b))

// creates row boundary
int row_line(int col)
{
    cout << endl;
    for (int i = 0; i < col; i++) {
        cout << " -----";
    }
    cout << endl;
}
 
// returns the count of alive neighbours
int count_live_neighbour_cell(const vector<vector<int> >& m, int r, int c)
{
    int i, j, count = 0;
    int row = m.size();
    int col = m[0].size();

    for (i = r - 1; i <= r + 1; i++) {
        for (j = c - 1; j <= c + 1; j++) {
            if ((i == r && j == c) || (i < 0 || j < 0)
                || (i >= row || j >= col)) {
                continue;
            }
            if (m[i][j] == 1) {
                count++;
            }
        }
    }
    return count;
}

// returns nextCanvas
vector<vector<int> > nextCanvas(const vector<vector<int> >& m)
{   
    int row = m.size();
    int col = m[0].size();

    vector<vector<int> > res (row, vector<int>(col));

    for (int r = 0; r < row; r++) {
        for (int c = 0; c < col; c++) {
            int neighbour_live_cell = count_live_neighbour_cell(m, r, c);

            if (m[r][c] == 1 && (neighbour_live_cell == 2 || neighbour_live_cell == 3))
                res[r][c] = 1;
            else if (m[r][c] == 0 && neighbour_live_cell == 3)
                res[r][c] = 1;
            else
                res[r][c] = 0;
        }
    }
    return res;
}

// read state matrix from file
void readMatrix(vector<vector<int> >& m){
    // Open the file
    std::ifstream inputFile("stateMatrix.txt");

    // Check if the file is opened successfully
    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file!" << std::endl;
        exit(1);
    }

    // Read the file line by line
    std::string line;
    while (std::getline(inputFile, line)) {
        // Vector to store the current row
        std::vector<int> row;

        // Process each character in the line
        for (size_t i = 0; i < line.length(); ++i) {
            // Convert character to integer and add to the row vector
            row.push_back(line[i]  - '0');
        }

        // Add the row vector to the matrix
        m.push_back(row);
    }

    // Close the file
    inputFile.close();
}

int main(int argc, char** argv)
{
    if(argc != 3){
        std::cerr << "Please specify your j and k in command line argument" << std::endl;
        std::cerr << "Usage: ./6-13 j k" << std::endl;
        exit(1);
    }

    int num_iter = atoi(argv[1]);
    int print_iter = atoi(argv[2]);
    int id, p;
    
    double elapsed_time;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier (MPI_COMM_WORLD);
    elapsed_time = - MPI_Wtime();

    // read matrix from file into state
    vector<vector<int> > matrix;
    readMatrix(matrix);
    
    int row = matrix.size();
    int col = matrix[0].size();

    // print matrix
    if(!id){
        cout << "Initial Stage:";
        row_line(col);
        for (int i = 0; i < row; i++) {
            cout << ":";
            for (int j = 0; j < col; j++) {
                cout << "  " << matrix[i][j] << "  :";
            }
            row_line(col);
        }
    }
    
    // calculate the index of submatrix of each processor
    int size;
    int smaller_size = row / p;
    int larger_size = smaller_size + 1;
    int num_larger_blocks = row % p;
    if (id < num_larger_blocks) size = larger_size;
    else size = smaller_size;
    int low_proc_value = (id*smaller_size + MIN(id, num_larger_blocks));
    int high_proc_value = low_proc_value + size - 1;

    // arguments for MPI_Gatherv
    std::vector<int> recvcounts(p);
    for (int i = 0; i < p; ++i) {
        if (i < num_larger_blocks) recvcounts[i] = larger_size * col;
        else recvcounts[i] = smaller_size * col;
    }
    std::vector<int> displs(p, 0);
    for (int i = 1; i < p; ++i) {
        displs[i] = displs[i - 1] + recvcounts[i - 1];
    }

    // portion of state in each processor
    vector<vector<int> > state;
    for (int i = low_proc_value; i < (low_proc_value + size); i++)
        state.push_back(matrix[i]);

    // start the iteration
    for(int ni = num_iter; ni > 0; ni -= print_iter){
        for(int pi = 0; pi < print_iter; ++pi){
            vector<int> receivedData1(col, 0), receivedData2(col, 0);
            // send neighbor processor bottom vector
            if(id != p-1){
                MPI_Send (state.back().data(), state.back().size(), MPI_INT, id+1, id, MPI_COMM_WORLD);
            }
            if(id){
                MPI_Recv (receivedData1.data(), receivedData1.size(), MPI_INT, id-1, id-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // send neighbor processor top vector
            if(id){
                MPI_Send (state.front().data(), state.front().size(), MPI_INT, id-1, p + id, MPI_COMM_WORLD);
            }
            if(id != p-1){
                MPI_Recv (receivedData2.data(), receivedData2.size(), MPI_INT, id+1, p + id+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            // concatenate these arrays
            if(id)
                state.insert(state.begin(), receivedData1); // need optimization
            if(id != p-1)
                state.push_back(receivedData2);

            // next canvas values based on live neighbour count
            state = nextCanvas(state);

            // restore the states
            if(id)
                state.erase(state.begin());
            if(id != p-1)
                state.pop_back();
        }

        // Flatten the 2D vector into a 1D array
        std::vector<int> sendBuffer(state.size() * state[0].size());
        int index = 0;
        for (int i = 0; i < state.size(); ++i) {
            for (int j = 0; j < state[0].size(); ++j) {
                sendBuffer[index++] = state[i][j];
            }
        }
        
        // gather the portion of states from all processors
        std::vector<int> recvBuffer(matrix.size() * matrix[0].size());
        MPI_Gatherv(sendBuffer.data(), sendBuffer.size(), MPI_INT, recvBuffer.data(), recvcounts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);

        // print next generation
        if(!id){
            cout << "\nNext Generation:";
            row_line(col);
            for (int i = 0; i < row; i++) {
                cout << ":";
                for (int j = 0; j < col; j++) {
                    cout << "  " << recvBuffer[i*col + j] << "  :";
                }
                row_line(col);
            }
        }
    }

    elapsed_time += MPI_Wtime();
    if(!id){
        std::cout << "elapsed_time: " << elapsed_time << std::endl;
    }

    MPI_Finalize();

    return 0;
}