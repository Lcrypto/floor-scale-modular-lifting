
/*
Copyright(c) 2014, Ilya Vorobyev und Vasiliy Usatyuk
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met :
*Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and / or other materials provided with the distribution.
* Neither the name of the <organization> nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include<iostream>
#include<string>
#include<algorithm>
#include<sstream>
#include".\myLib\localOptimization.h"


using namespace std;
#define pii pair<int, int>
#define ll long long


int girthFirst = 0;
int aceFirst = 0;

string toStr(ll x) {
    stringstream ss;
    ss << x;
    return ss.str();
}

bool checkInput(string MTR, string SIZES) {
    bool validInput = 1;
    if (MTR == "") {
        validInput = 0;
        cerr << "wrong matrix file\n";
        cerr << "example: -matrix matrix.txt\n";
        cerr << "matrix format:\nvar\tcheck\tcirc\na_{1,1} a_{1,2}\t...\ta_{1,var}\n...\n...\n...\na_{check,1} a_{check,2}\t...\ta_{check,var}\n";
        cerr << endl;
    }
    if (SIZES == "") {
        validInput = 0;
        cerr << "wrong sizes\n";
        cerr << "example: -sizes sizes.txt\n";
        cerr << "sizes format:\nvar1\tcheck1\tgirth1\tit1\n.\n.\n.\nvar_n\tcheck_n\tgirth_n\tit_n";
        cerr << endl;
    }
    return validInput;
}

bool nextCombination(vector<int>& a, int n) {
    int k = (int)a.size();
    for (int i = k - 1; i >= 0; --i)
        if (a[i] < n - k + i + 1) {
            ++a[i];
            for (int j = i + 1; j < k; ++j)
                a[j] = a[j - 1] + 1;
            return true;
        }
    return false;
}



//馬耶夫斯基同性戀
//for regular
vector<int> getPermanent(const vector<vector<int> >& a, vector<int>& used, int mod, int step = 0) {
    int n = a.size();
    vector<int> res(mod, 0);
    if (step + 1 == n) {
        for (int i = 0; i < n + 1; ++i) {
            if (used[i])
                continue;
            if (a[n - 1][i] == -1)
                return res;
            res[a[n - 1][i]] = 1;
            return res;
        }
    }
    for (int i = 0; i < n + 1; ++i) {
        if (used[i])
            continue;
        if (a[step][i] == -1)
            continue;
        used[i] = 1;
        vector<int> cur = getPermanent(a, used, mod, step + 1);
        for (int j = 0; j < cur.size(); ++j) {
            int x = j + a[step][i];
            if (x >= mod)
                x -= mod;
            res[x] ^= cur[j];
        }
        used[i] = 0;
    }
    return res;
}



//馬耶夫斯基同性戀
//for irregular
vector<int> getPermanent(const vector<vector<vector<int> > >& a, vector<int>& used, int mod, int step = 0) {
    int n = a.size();
    vector<int> res(mod, 0);
    if (step + 1 == n) {
        for (int i = 0; i < n + 1; ++i) {
            if (used[i])
                continue;
            for (int j = 0; j < a[n - 1][i].size(); ++j)
                res[a[n - 1][i][j]] = 1;
            return res;
        }
    }
    for (int i = 0; i < n + 1; ++i) {
        if (used[i])
            continue;
        if (a[step][i].empty())
            continue;
        used[i] = 1;
        vector<int> cur = getPermanent(a, used, mod, step + 1);
        for (int j = 0; j < cur.size(); ++j) {
            if (cur[j] == 0)
                continue;
            for (int jj = 0; jj < a[step][i].size(); ++jj) {
                int x = j + a[step][i][jj];
                if (x >= mod)
                    x -= mod;
                res[x] ^= cur[j];
            }
        }
        used[i] = 0;
    }
    return res;
}

//for regular
int getWeight(const vector<vector<int> >& a, int erInd, int mod) {
    int n = a.size();
    vector<int> used(n + 1, 0);
    used[erInd] = 1;
    vector<int> pol = getPermanent(a, used, mod);
    int res = 0;
    for (int i = 0; i < pol.size(); ++i)
        res += pol[i];
    return res;

}

//for irregular
int getWeight(const vector<vector<vector<int> > >& a, int erInd, int mod) {
    int n = a.size();
    vector<int> used(n + 1, 0);
    used[erInd] = 1;
    vector<int> pol = getPermanent(a, used, mod);
    int res = 0;
    for (int i = 0; i < pol.size(); ++i)
        res += pol[i];
    return res;
}

//for regular
int solve(const vector<int>& mask, const vector<vector<int> >& mtr, int mod) {
    vector<vector<int> > newMtr(mtr.size(), vector<int>(mask.size()));
    for (int i = 0; i < mtr.size(); ++i)
        for (int j = 0; j < mask.size(); ++j) {
            newMtr[i][j] = mtr[i][mask[j]];
        }
    int res = 0;
    for (int i = 0; i < mask.size(); ++i)
        res += getWeight(newMtr, i, mod);
    return res;
}

//for irregular
int solve(const vector<int>& mask, const vector<vector<vector<int> > >& mtr, int mod) {
    vector<vector<vector<int> > > newMtr(mtr.size(), vector<vector<int> >(mask.size()));
    for (int i = 0; i < mtr.size(); ++i)
        for (int j = 0; j < mask.size(); ++j) {
            newMtr[i][j] = mtr[i][mask[j]];
        }
    int res = 0;
    for (int i = 0; i < mask.size(); ++i)
        res += getWeight(newMtr, i, mod);
    return res;
}

//for regular
int countBound(const vector<vector<int> > & mtr, int mod) {
    int J = mtr.size(), I = mtr[0].size();
    if (I <= J)
        return -1;
    vector<int> mask(J + 1, 0);
    for (int i = 0; i < J + 1; ++i)
        mask[i] = i;
    int res = -1;
    do {
        int cur = solve(mask, mtr, mod);
        if (cur > 0) {
            if ((res < 0) || (cur < res))
                res = cur;
        }
    } while (nextCombination(mask, I - 1));
    return res;
}

//for irregular
int countBound(const vector<vector<vector<int> > > & mtr, int mod) {
    int J = mtr.size(), I = mtr[0].size();
    if (I <= J)
        return -1;
    vector<int> mask(J + 1, 0);
    for (int i = 0; i < J + 1; ++i)
        mask[i] = i;
    int res = -1;
    do {
        int cur = solve(mask, mtr, mod);
        if (cur > 0) {
            if ((res < 0) || (cur < res))
                res = cur;
        }
    } while (nextCombination(mask, I - 1));
    return res;
}

void print(const vector<int>& a) {
    for (int i = 0; i < a.size(); ++i)
        cout << a[i] << " ";
    cout << endl;
}


vector<int> parse(string s) {
    if ((s.empty()) || (s[0] == '-'))
        return vector<int>();
    vector<int> res;
    int x = 0;
    for (int i = 0; i < s.size(); ++i) {
        if (s[i] == '&') {
            res.push_back(x);
            x = 0;
        }
        else
            x = 10 * x + (s[i] - '0');
    }
    res.push_back(x);
    return res;
}

vector<vector<int> > getSmallProto(const vector<vector<int> >& proto, int var, int check) {
    vector<vector<int> > res(check, vector<int>(var));
    for (int i = 0; i < check; ++i)
        for (int j = 0; j < var; ++j)
            res[i][j] = proto[i][j];
    return res;
}


pair<int, int> getGirthAndAce(const vector<vector<vector<int> > >& a, int circ, const vector<vector<int> >& proto, int upGirth) {
    vector<int> degVar(a[0].size(), 0);

    for (int r = 0; r < a.size(); ++r) {
        for (int c = 0; c < a[r].size(); ++c) {
            degVar[c] += a[r][c].size();
        }
    }

    pair<int, int> girthAce(-1, -1);
    for (int g = 4; g <= upGirth; g += 2) {
        CycleEnum enumerator(g, proto);
        if (!enumerator.init()) {
            continue;
        }
        do {
            vector<entry> cycle = enumerator.cycle;
            int res = 0;
            for (int i = 0; i < g; ++i) {
                int x = a[cycle[i].r][cycle[i].c][cycle[i].id];
                if (i & 1)
                    res += x;
                else
                    res -= x;
            }
            res = (((res % circ) + circ) % circ);
            if (res == 0) {
                if (girthAce.first == -1)
                    girthAce.first = g;
                int curAce = 0;
                for (int i = 0; i < cycle.size(); ++i) {
                    if (i & 1)
                        curAce += degVar[cycle[i].c] - 2;
                }
                if ((girthAce.second == -1) || (curAce < girthAce.second))
                    girthAce.second = curAce;
            }
        } while (enumerator.next());
    }
    if (girthAce.first == -1)
        girthAce.first = 100;
    if (girthAce.second == -1)
        girthAce.second = 1000000;
    return girthAce;
}


struct record {
    int scaleFactor, girth, ace, dist;
    bool operator<=(const record& rec) {
        if (girthFirst) {
            return (girth < rec.girth) || ((girth == rec.girth) && (ace <= rec.ace));
        }
        if (aceFirst) {
            return (ace < rec.ace) || ((ace == rec.ace) && (girth <= rec.girth));
        }
        return (girth <= rec.girth) && (ace <= rec.ace) && (dist <= rec.dist);
    }
};

ostream& operator<<(ostream& os, const record& rec) {
    os << "scale factor = " << rec.scaleFactor << "\tgirth = " << rec.girth << "\tace = " << rec.ace;
    if (rec.dist != 1000) {
        os << "\tdist from th7 = " << rec.dist;
    }
    os << endl;
    return os;
}

int main(int argc, char* argv[]) {
    //Initialization
    bool validInput = 1;
    vector<vector<int> > PROTOGRAPH;
    int CIRCULANT_SIZE = -1;
    ll VARIABLE_NODES;
    ll CHECK_NODES;
    string MTR_FILE = "";
    string SIZES_FILE = "";
    bool th7 = 0;
    int distUp = 1000;
    int upGirth = 10;
    
    //reading parameters from cmd
    for (int i = 1; i < argc; ++i) {
        if (string(argv[i]) == "-matrix") {
            MTR_FILE = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-sizes") {
            SIZES_FILE = argv[i + 1];
            ++i;
            continue;
        }
        if (string(argv[i]) == "-th7") {
            th7 = 1;
            continue;
        }
        if (string(argv[i]) == "-distUp") {
            string strDist = argv[i + 1];
            ++i;
            stringstream sstrDist(strDist);
            sstrDist >> distUp;
            continue;
        }
        if (string(argv[i]) == "-upGirth") {
            string strGirth = argv[i + 1];
            ++i;
            stringstream sstrGirth(strGirth);
            sstrGirth >> upGirth;
            continue;
        }
        if (string(argv[i]) == "-girthFirst") {
            girthFirst = 1;
            continue;
        }
        if (string(argv[i]) == "-aceFirst") {
            aceFirst = 1;
            continue;
        }
    }
    if ((!checkInput(MTR_FILE, SIZES_FILE)) || (!validInput))
        return 0;

    //reading matrix -- START
    freopen(MTR_FILE.c_str(), "r", stdin);
    cin >> VARIABLE_NODES >> CHECK_NODES >> CIRCULANT_SIZE;
    vector<vector<vector<int> > > mtr(CHECK_NODES, vector<vector<int> >(VARIABLE_NODES));
    PROTOGRAPH.assign(CHECK_NODES, vector<int>(VARIABLE_NODES, 0));
    for (int i = 0; i < CHECK_NODES; ++i) {
        for (int j = 0; j < VARIABLE_NODES; ++j) {
            string toParse;
            cin >> toParse;
            mtr[i][j] = parse(toParse);
            PROTOGRAPH[i][j] = mtr[i][j].size();
        }
    }
    fclose(stdin);


    //reading matrix -- END

    freopen(SIZES_FILE.c_str(), "r", stdin);
    int numberOfSizes;
    cin >> numberOfSizes;
    vector<int> varNodes(numberOfSizes), checkNodes(numberOfSizes), circ(numberOfSizes), numOfIter(numberOfSizes);// , resGirth(numberOfSizes), resACE(numberOfSizes), scaling(numberOfSizes);
    vector<vector<record> > records(numberOfSizes);
    for (int i = 0; i < numberOfSizes; ++i) {
        cin >> varNodes[i] >> checkNodes[i] >> circ[i] >> numOfIter[i];
    }
    fclose(stdin);

    //folder and filenames generation
    string folderName = "scaling" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE);
    string outputFilename = folderName + "/" + toStr(VARIABLE_NODES) + "_" + toStr(CHECK_NODES) + "_" + toStr(CIRCULANT_SIZE);
    outputFilename += "matrix_from_" + MTR_FILE;
    outputFilename += "sizes_from_" + SIZES_FILE;
    //system(("mkdir " + folderName).c_str());

    for (int sizeId = 0; sizeId < numberOfSizes; ++sizeId) {
        cerr << "start to compute scale factor " << sizeId + 1 << "(circulat = " << circ[sizeId] << ")" << endl;
        vector<vector<vector<int> > > newMtr(checkNodes[sizeId], vector<vector<int> >(varNodes[sizeId]));
        for (int rIt = 1; rIt <= min(numOfIter[sizeId], CIRCULANT_SIZE); ++rIt) {
            //cerr << rIt << endl;
            int r = rIt;
            if (numOfIter[sizeId] < CIRCULANT_SIZE)
                r = rand() % CIRCULANT_SIZE;
            record cur;
            cur.scaleFactor = r;
            cur.dist = distUp;
            bool okScale = 1;
            for (int i = 0; i < checkNodes[sizeId]; ++i) {
                for (int j = 0; j < varNodes[sizeId]; ++j) {
                    newMtr[i][j].resize(PROTOGRAPH[i][j]);
                    for (int nodeId = 0; nodeId < mtr[i][j].size(); ++nodeId) {

                        newMtr[i][j][nodeId] = (((mtr[i][j][nodeId] * r) % CIRCULANT_SIZE) * circ[sizeId]) / CIRCULANT_SIZE;
                    }
                    
                    for (int i1 = 0; i1 < mtr[i][j].size(); ++i1) {
                        for (int i2 = 0; i2 < i1; ++i2) {
                            if (newMtr[i][j][i1] == newMtr[i][j][i2]) {
                                okScale = 0;
                                break;
                            }
                        }
                    }
                }
            }
            if (!okScale)
                continue;


            pair<int, int> girthAce = getGirthAndAce(newMtr, circ[sizeId], getSmallProto(PROTOGRAPH, varNodes[sizeId], checkNodes[sizeId]), upGirth);
            cur.girth = girthAce.first;
            cur.ace = girthAce.second;
            bool goodRec = 1;
            for (int i = 0; i < records[sizeId].size(); ++i) {
                if (cur <= records[sizeId][i]) {
                    goodRec = 0;
                    break;
                }
            }
            if (!goodRec)
                continue;
            if (th7) {
                cur.dist = countBound(newMtr, circ[sizeId]);
            }
            for (int i = 0; i < records[sizeId].size(); ++i) {
                if (cur <= records[sizeId][i]) {
                    goodRec = 0;
                    break;
                }
            }
            if (goodRec) {
                records[sizeId].push_back(cur);
                cerr << "new record for matrix " << sizeId + 1 << ":\t" << cur;
            }

        }
   
    }
    //freopen("out.txt", "w", stdout);
    FILE* out = fopen("detailedOutput.txt", "w");
    for (int i = 0; i < records.size(); ++i) {
        fprintf(out, "var = %d\tcheck = %d\tcirc = %d\n", varNodes[i], checkNodes[i], circ[i]);
        cerr << "var = " << varNodes[i] << "\tcheck = " << checkNodes[i] << "\tcirc = " << circ[i] << endl;
        for (int j = 0; j < records[i].size(); ++j) {
            bool goodRec = 1;
            for (int k = j + 1; k < records[i].size(); ++k) {
                if (records[i][j] <= records[i][k])
                    goodRec = 0;
            }
            if (!goodRec)
                continue;
            cerr << records[i][j];
            if (th7)
                fprintf(out, "scale factor = %d\tgirth = %d\tace on cycles of length <= %d = %d\tdist from th7 = %d\n", records[i][j].scaleFactor, records[i][j].girth, upGirth, records[i][j].ace, records[i][j].dist);
            else
                fprintf(out, "scale factor = %d\tgirth = %d\tace on cycles of length <= %d = %d\n", records[i][j].scaleFactor, records[i][j].girth, upGirth, records[i][j].ace);
        }
    }
    fclose(out);

    if (!th7) {
        out = fopen("coef.txt", "w");
        fprintf(out, "%d\n", records.size());
        for (int i = 0; i < records.size(); ++i) {
            if (records[i].empty()) {
                cerr << "there is no good scale factor\n";
            }
            
            for (int j = 0; j < records[i].size(); ++j) {
                bool goodRec = 1;
                for (int k = j + 1; k < records[i].size(); ++k) {
                    if (records[i][j] <= records[i][k])
                        goodRec = 0;
                }
                if (!goodRec)
                    continue;
                fprintf(out, "%d\t%d\n", circ[i], records[i][j].scaleFactor);
                
            }
        

        }
        fclose(out);

    }
    else {
        out = fopen("coef.txt", "w");
        for (int i = 0; i < records.size(); ++i) {
            bool first = 1;
            for (int j = 0; j < records[i].size(); ++j) {
                bool goodRec = 1;
                for (int k = j + 1; k < records[i].size(); ++k) {
                    if (records[i][j] <= records[i][k])
                        goodRec = 0;
                }
                if (!goodRec)
                    continue;
                if (first)
                    fprintf(out, "%d", records[i][j].scaleFactor);
                else
                    fprintf(out, ",%d", records[i][j].scaleFactor);
                first = 0;
            }
            fprintf(out, "\n");

        }
        fclose(out);
    }
    return 0;
}
