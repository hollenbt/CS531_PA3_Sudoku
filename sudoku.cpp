#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;

class Inferer {
    int **board;

public:

    // With at most 9 elements in each set, the asymptotic time complexity advantage
    // of unordered_set vs. set is irrelevant. We can experiment with both to see
    // which has better time and space performance.

    // Cells values: if assigned, set will be emtpy
    vector<vector<set<int>>> cell;
    // Value locations: if assigned, set will be empty
    vector<vector<set<int>>> row;              // row[row #][value]:     set of col #'s
    vector<vector<set<int>>> col;              // col[col #][value]:     set of row #'s
    vector<vector<set<pair<int,int>>>> block; // block[block #][value]: set of coords


    Inferer(int **b) : board(b) {}

    void print() {
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j)
                cout << ' ' << board[i][j];
            cout << endl;
        }
        cout << endl;
    }

    bool initialize() {
        reset();
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                int v = board[i][j];
                if (0 < v && v <= 9) {
                    if (assign(i, j, v) == 2) return false;
                }
            }
        }
        return true;
    }

    void reset() {
        cell  = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 1 ,2, 3, 4, 5, 6, 7, 8, 9 })));
        row   = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 0, 1 ,2, 3, 4, 5, 6, 7, 8 })));
        col   = vector<vector<set<int>>>(9, vector<set<int>>(9, set<int>({ 0, 1 ,2, 3, 4, 5, 6, 7, 8 })));
        block = vector<vector<set<pair<int,int>>>>(9, vector<set<pair<int, int>>>(9, set<pair<int, int>>()));
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) {
                pair<int, int> p(i, j);
                int b = get_block(p);
                for (int k = 0; k < 9; ++k)
                    block[b][k].insert(p);
            }

    }

    int get_block(const pair<int, int> &p) {
        return 3 * (p.first / 3) + (p.second / 3);
    }

    /***********************
     * Variable Assignment *
     ***********************/
    // Return values:
    //   1 - successful assignment
    //   2 - conflict (i.e. empty domain for unassigned variable), need to backtrack
    int assign(int i, int j, int v) {
        // Check whether assignment is valid
        if (cell[i][j].find(v) == cell[i][j].end())
            return 2;

        pair<int, int> p(i, j);
        int b = get_block(p);

        // Clear the sets to indicate assignment
        // (and prevent domain restriction from
        // returning empty domain error)
        cell[i][j].clear();
        row[i][v - 1].clear();
        col[j][v - 1].clear();
        block[b][v - 1].clear();

        // Remove value from rest of row and column.
        // Remove coordinates from all other values.
        for (int k = 0; k < 9; ++k) {
            if (eliminate(i, k, v)) return 2;
            if (eliminate(k, j, v)) return 2;
            if (eliminate(i, j, k + 1)) return 2;
        }

        // Remove value from rest of block
        int bi = 3 * (i / 3);
        int bj = 3 * (j / 3);
        for (int x = bi; x < bi + 3; ++x)
            for (int y = bj; y < bj + 3; ++y)
                if(eliminate(x, y, v)) return 2;

        // Make assignment on board
        board[i][j] = v;
        // could print board here

        return 1;
    }

    /**********************
     * Domain restriction *
     **********************/
    // Return: true if restriction makes a domain empty,
    //         false otherwise
    bool eliminate(int i, int j, int v) {
        if (cell[i][j].erase(v) && cell[i][j].empty()) return true;
        if (row[i][v - 1].erase(j) && row[i][v - 1].empty()) return true;
        if (col[j][v - 1].erase(i) && col[j][v - 1].empty()) return true;
        pair<int, int> p(i, j);
        int b = get_block(p);
        if (block[b][v - 1].erase(p) && block[b][v - 1].empty()) return true;
        return false;
    }

    /*******************
     * Inference Rules *
     *******************/
    // Return values:
    //   0 - no application found
    //   1 - application found and applied
    //   2 - conflict, need to backtrack

    int naked_single() {
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                if (cell[i][j].size() == 1) {
                    cout << "Naked single found in cell (" << i << ", " << j << ')' << endl;
                    return assign(i, j, *cell[i][j].begin());
                }
        return 0;
    }

    int hidden_single() {
        for (int x = 0; x < 9; ++x) {
            for (int v = 0; v < 9; ++v) {
                if (row[x][v].size() == 1) {
                    cout << "Hidden single found in cell (" << x << ", " << *row[x][v].begin() << ')' << endl;
                    return assign(x, *row[x][v].begin(), v + 1);
                }
                if (col[x][v].size() == 1) {
                    cout << "Hidden single found in cell (" << *col[x][v].begin() << ", " << x << ')' << endl;
                    return assign(*col[x][v].begin(), x, v + 1);
                }
                if (block[x][v].size() == 1) {
                    pair<int, int> p = *block[x][v].begin();
                    cout << "Hidden single found in cell (" << p.first << ", " << p.second << ')' << endl;
                    return assign(p.first, p.second, v + 1);
                }
            }
        }
        return 0;
    }

    // If we are performing this, there are no naked/hidden singles, so all row/col/block.size() > 1
    // (if referenced by a cell, otherwise could be 0). This means we only need to check that the
    // larger cell set size is > 2, and not that the smaller cell set size is >= 2.
    int naked_pair() {
        // Check by row
        for (int r = 0; r < 9; ++r) {
            for (int c1 = 0; c1 < 8; ++c1) {
                // Check for cells that have only 2 values
                if (cell[r][c1].size() == 2) {
                    set<int>::iterator it = cell[r][c1].begin();
                    bool second = false;
                    int n1 = row[r][*it - 1].size();
                    ++it;
                    int n2 = row[r][*it - 1].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the values must appear in more than 2 cells in row
                        it = cell[r][c1].begin();
                        if (second) ++it; // Iterate through smaller of the two column sets
                        for (set<int>::iterator cit = row[r][*it - 1].upper_bound(c1); cit != row[r][*it - 1].end(); ++cit)
                            if (cell[r][c1] == cell[r][*cit]) { // If naked pair found
                                cout << "Row: Naked pair found in cells (" << r << ", " << c1 << ')';
                                cout << " and (" << r << ", " << *cit << ")" << endl;
                                // Eliminate the pair of values from all other cells in the row
                                for (int v : cell[r][c1]) {
                                    for (int c : row[r][v - 1])
                                        if (c != c1 && c != *cit) {
                                            cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                            if (eliminate(r, c, v)) return 2;
                                        }
                                }
                                cout << endl;
                                return 1;
                            }
                    }
                }
            }
        }

        // Check by column
        for (int c = 0; c < 9; ++c) {
            for (int r1 = 0; r1 < 8; ++r1) {
                // Check for cells that have only 2 values
                if (cell[r1][c].size() == 2) {
                    set<int>::iterator it = cell[r1][c].begin();
                    bool second = false;
                    int n1 = col[c][*it - 1].size();
                    ++it;
                    int n2 = col[c][*it - 1].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the values must appear in more than 2 cells in column
                        it = cell[r1][c].begin();
                        if (second) ++it; // Iterate through smaller of the two row sets
                        for (set<int>::iterator rit = col[c][*it - 1].upper_bound(r1); rit != col[c][*it - 1].end(); ++rit)
                            if (cell[r1][c] == cell[*rit][c]) { // If naked pair found
                                cout << "Column: Naked pair found in cells (" << r1 << ", " << c << ')';
                                cout << " and (" << *rit << ", " << c << ")" << endl;
                                // Eliminate the pair of values
                                for (int v : cell[r1][c]) {
                                    // From all other cells in the column
                                    for (int r : col[c][v - 1])
                                        if (r != r1 && r != *rit) {
                                            cout << "Eliminating (" << r << ", " << c << ", " << v << ")" << endl;
                                            if (eliminate(r, c, v)) return 2;
                                        }
                                }
                                cout << endl;
                                return 1;
                            }
                    }
                }
            }
        }

        // Check by block
        for (int b = 0; b < 9; ++b) {
            int r_beg = 3 * (b / 3);
            int c_beg = 3 * (b % 3);
            for (int r1 = r_beg; r1 < r_beg + 3; ++r1) {
                for (int c1 = c_beg; c1 < c_beg + 3; ++c1) {
                    // Check for cells that have only 2 values
                    if (cell[r1][c1].size() == 2) {
                        pair<int, int> p1(r1, c1);
                        set<int>::iterator it = cell[r1][c1].begin();
                        bool second = false;
                        int n1 = block[b][*it - 1].size();
                        ++it;
                        int n2 = block[b][*it - 1].size();
                        if (n1 > n2) {
                            swap(n1, n2);
                            second = true;
                        }
                        if (n2 > 2) { // One of the values must appear in more than 2 cells in block
                            it = cell[r1][c1].begin();
                            if (second) ++it; // Iterate through smaller of the two coordinate sets
                            for (set<pair<int, int>>::iterator pit = block[b][*it - 1].upper_bound(p1); pit != block[b][*it - 1].end(); ++pit)
                                if (cell[r1][c1] == cell[pit->first][pit->second]) { // If naked pair found
                                    cout << "Block: Naked pair found in cells (" << r1 << ", " << c1 << ')';
                                    cout << " and (" << pit->first << ", " << pit->second << ")" << endl;
                                    // Eliminate the pair of values from all other cells in the block
                                    for (int v : cell[r1][c1]) {
                                        for (const pair<int, int> &p : block[b][v - 1])
                                            if (p != p1 && p != *pit) {
                                                cout << "Eliminating (" << p.first << ", " << p.second << ", " << v << ")" << endl;
                                                if (eliminate(p.first, p.second, v)) return 2;
                                            }
                                    }
                                    cout << endl;
                                    return 1;
                                }
                        }
                    }
                }
            }
        }

        return 0;
    }

    // If we are performing this, there are no naked/hidden singles or naked pairs, so
    // we can ignore pairs of cells where the larger value set size is 2 or less.
    int hidden_pair() {
        for (int x = 0; x < 9; ++x) {
            for (int v = 0; v < 8; ++v) {

                // Check row for values that appear in only 2 cells
                if (row[x][v].size() == 2) {
                    set<int>::iterator it = row[x][v].begin();
                    bool second = false;
                    int n1 = cell[x][*it].size();
                    ++it;
                    int n2 = cell[x][*it].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values
                        it = row[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[x][*it].upper_bound(v + 1); vit != cell[x][*it].end(); ++vit)
                            if (row[x][v] == row[x][*vit - 1]) { // If hidden pair found
                                cout << "Row: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                // Eliminate all other values from the pair of cells.
                                for (int c : row[x][v]) {
                                    for (int val : cell[x][c])
                                        if (val != v + 1 && val != *vit) {
                                            cout << "Eliminating (" << x << ", " << c << ", " << val << ")" << endl;
                                            if (eliminate(x, c, val)) return 2;
                                        }
                                }
                                return 1;
                            }
                    }
                }

                // Check column for values that appear in only 2 cells
                if (col[x][v].size() == 2) {
                    set<int>::iterator it = col[x][v].begin();
                    bool second = false;
                    int n1 = cell[*it][x].size();
                    ++it;
                    int n2 = cell[*it][x].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values
                        it = col[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[*it][x].upper_bound(v + 1); vit != cell[*it][x].end(); ++vit)
                            if (col[x][v] == col[x][*vit - 1]) { // If hidden pair found
                                cout << "Column: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                // Eliminate all other values from the pair of cells.
                                for (int r : col[x][v]) {
                                    for (int val : cell[r][x])
                                        if (val != v + 1 && val != *vit) {
                                            cout << "Eliminating (" << r << ", " << x << ", " << val << ")" << endl;
                                            if (eliminate(r, x, val)) return 2;
                                        }
                                }
                                return 1;
                            }
                    }
                }

                // Check block for values that appear in only 2 cells
                if (block[x][v].size() == 2) {
                    set<pair<int, int>>::iterator it = block[x][v].begin();
                    bool second = false;
                    int n1 = cell[it->first][it->second].size();
                    ++it;
                    int n2 = cell[it->first][it->second].size();
                    if (n1 > n2) {
                        swap(n1, n2);
                        second = true;
                    }
                    if (n2 > 2) { // One of the cells must have more than 2 values
                        it = block[x][v].begin();
                        if (second) ++it; // Iterate through smaller of the two value sets
                        for (set<int>::iterator vit = cell[it->first][it->second].upper_bound(v + 1); vit != cell[it->first][it->second].end(); ++vit)
                            if (block[x][v] == block[x][*vit - 1]) { // If hidden pair found
                                cout << "Block: Hidden pair found with values " << v + 1 << ", " << *vit << endl;
                                // Eliminate all other values from the pair of cells.
                                for (const pair<int, int> &p : block[x][v]) {
                                    for (int val : cell[p.first][p.second])
                                        if (val != v + 1 && val != *vit) {
                                            cout << "Eliminating (" << p.first << ", " << p.second << ", " << val << ")" << endl;
                                            if (eliminate(p.first, p.second, val)) return 2;
                                        }
                                }
                                return 1;
                            }
                    }
                }
            }
        }
        return 0;
    }

    int naked_triple() {
        return 0;
    }

    int hidden_triple() {
        return 0;
    }

    // Inference Schemes
    int no_inference() {
        return 0;
    }

    int up_to_singles() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0)
                r = hidden_single();
        }
        return r;
    }

    int up_to_pairs() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0) {
                r = hidden_single();
                if (r == 0) {
                    r = naked_pair();
                    if (r == 0)
                        r = hidden_pair();
                }
            }
        }
        return r;
    }

    int up_to_triples() {
        int r = 1;
        while (r == 1) {
            r = naked_single();
            if (r == 0) {
                r = hidden_single();
                if (r == 0) {
                    r = naked_pair();
                    if (r == 0) {
                        r = hidden_pair();
                        if (r == 0) {
                            r = naked_triple();
                            if (r == 0)
                                r == hidden_triple();
                        }
                    }
                }
            }
        }
        return r;
    }
};


// Generate row, column, block, and cell sets
// Apply inference rules repeatedly until none apply
// 


// choose next square with fixed baseline

// choose next square with most restricted variable


// backtracking search
void backtracking_search() {

}

int main() {
    int **board = new int*[9];
    for (int i = 0; i < 9; ++i)
        board[i] = new int[9];

    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            cin >> board[i][j];

    Inferer inferer(board);
    inferer.print();
    inferer.initialize();
    inferer.up_to_pairs();
    inferer.print();

    for (int i = 0; i < 9; ++i)
        delete[] board[i];
    delete[] board;

    return 0;
}