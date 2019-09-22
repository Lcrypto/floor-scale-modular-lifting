#pragma once
#include"..\myLib\CycleEnum.h"
#include"irregularLDPC.h"

//瑪雅混蛋
struct Tiii {
    int first, second, third;
    Tiii() {}
    Tiii(int _first, int _second, int _third) : first(_first), second(_second), third(_third) {}
    Tiii(const entry& e) : first(e.r), second(e.c), third(e.id) {}
    bool operator==(const Tiii& rhs) const {
        return (first == rhs.first) && (second == rhs.second) && (third == rhs.third);
    }
    bool operator!=(const Tiii& rhs) const {
        return (first != rhs.first) || (second != rhs.second) || (third != rhs.third);
    }
    bool operator<(const Tiii& rhs) const {
        return (first < rhs.first) ||
            ((first == rhs.first) && (second < rhs.second)) ||
            ((first == rhs.first) && (second == rhs.second) && (third == rhs.third));
    }
};

struct submatrix {
    int var, check, circ, girth, ace;
};

class LocalOptParallel {
    struct Cycle {
        vector<Tiii> cycle;
        int uniqueNodes = 1;
        bool oneEntry = true;
        Cycle(const vector<entry>& _cycle) {
            cycle.resize(_cycle.size());
            for (size_t i = 0; i < cycle.size(); ++i)
                cycle[i] = _cycle[i];
            for (int i = 1; i < cycle.size(); ++i) {
                if (cycle[i] != cycle[i - 1]) {
                    oneEntry = false;
                    break;
                }
            }
            for (int i = 1; i < cycle.size(); ++i) {
                ++uniqueNodes;
                for (int j = 0; j < i; ++j)
                    if (cycle[i] == cycle[j]) {
                        --uniqueNodes;
                        break;
                    }
            }
        }
    };
    int checkNodes, variableNodes;
    vector<vector<int> > varDeg, checkDeg;
    vector<vector<int> > protograph;
    vector<vector<vector<int> > > submatrixProto;
    vector<submatrix> target;
    vector<vector<int> > maxGirth;
    int maxCirc = 0;
    vector<vector<vector<int> > > mtr;
    bool noCycles = false;
    vector<int> scalingFactor;
    ll totalNumberOfCycles;

public:

    LocalOptParallel(vector<submatrix> _target, const vector<vector<vector<int> > >& _mtr) {
        checkNodes = _mtr.size();
        variableNodes = _mtr[0].size();
        mtr.assign(checkNodes, vector<vector<int> >(variableNodes));
        protograph.assign(checkNodes, vector<int>(variableNodes));
        mtr = _mtr;
        target = _target;
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                protograph[r][c] = mtr[r][c].size();
            }
        }

        submatrixProto.resize(target.size());
        varDeg.assign(target.size(), vector<int>(variableNodes, 0));
        checkDeg.assign(target.size(), vector<int>(checkNodes, 0));

        for (int targetId = 0; targetId < target.size(); ++targetId) {
            submatrixProto[targetId].resize(target[targetId].check);
            for (int r = 0; r < target[targetId].check; ++r) {
                submatrixProto[targetId][r].resize(target[targetId].var);
                for (int c = 0; c < target[targetId].var; ++c) {
                    submatrixProto[targetId][r][c] = protograph[r][c];
                    varDeg[targetId][c] += protograph[r][c];
                    checkDeg[targetId][r] += protograph[r][c];

                }
            }
        }

        maxGirth.assign(checkNodes, vector<int>(variableNodes, 0));
        for (int i = 0; i < target.size(); ++i) {
            int curGirth = target[i].girth - 2;
            if (target[i].ace > 0)
                curGirth = target[i].girth;
            for (int r = 0; r < target[i].check; ++r) {
                for (int c = 0; c < target[i].var; ++c) {
                    maxGirth[r][c] = max(maxGirth[r][c], curGirth);
                }
            }
            maxCirc = max(maxCirc, target[i].circ);
        }
        scalingFactor.assign(target.size(), 1);
    }

    //returns number of remaining cycles
    ll optimize(bool scaling = 0, bool coupling = 0) {
        computeTotalNumberOfCycles(coupling);
        long long numberOfValuesToAssign = 0;
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                numberOfValuesToAssign += protograph[r][c];
            }
        }
        int movesWithoutChange = 0;
        bool cycleExist = 0;
        long long moves = 0;
        //vector<ll> numCycMem;
        ll numberOfActCycles = getTotalNumberOfActiveCycles(mtr, coupling);
        ll minNumOfActiveCycles = numberOfActCycles;

        while (true) {
            for (int r = 0; r < checkNodes; ++r) {
                for (int c = 0; c < variableNodes; ++c) {
                    for (int id = 0; id < protograph[r][c]; ++id) {
                        //cerr << r << " " << c << endl;
                        //int r = getRand(checkNodes);
                        //int c = getRand(variableNodes);
                        if (protograph[r][c] == 0)
                            continue;
                        ll curNumOfActiveCycles = numberOfActCycles;//getTotalNumberOfActiveCycles(mtr, coupling);

                        ++movesWithoutChange;
                        ++moves;
                        if (numberOfActCycles == 0)
                            return 0;
                        if (movesWithoutChange > 1 + numberOfValuesToAssign) {
                           
                            return curNumOfActiveCycles;

                        }
                        vector<ll> numberOfCycles(maxCirc, 0);
                        for (int targetId = 0; targetId < target.size(); ++targetId) {
                            if (r >= target[targetId].check)
                                continue;
                            if (c >= target[targetId].var)
                                continue;
                            int curGirth = target[targetId].girth - 2;
                            if (target[targetId].ace > 0)
                                curGirth += 2;
                            for (int girth = 4; girth <= curGirth; girth += 2) {
                                //cerr << r << " " << c << " " << girth << endl;
                                CycleEnum enumerator(girth, submatrixProto[targetId]);
                                if (!enumerator.init(r, c, id))
                                    continue;
                                do {
                                    Cycle cycle(enumerator.cycle);
                                    if ((enumerator.cycle[0].r != r) || (enumerator.cycle[0].c != c))
                                        break;
                                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                        continue;
                                    processCycle(cycle, numberOfCycles, mtr, target[targetId], scalingFactor[targetId]);
                                } while (enumerator.next());
                            }
                        }
                        for (int i = 0; i < maxCirc; ++i) {
                            if (i == mtr[r][c][id])
                                continue;
                            if (numberOfCycles[mtr[r][c][id]] == 0)
                                break;
                            if (numberOfCycles[i] < numberOfCycles[mtr[r][c][id]]) {
                                numberOfActCycles += numberOfCycles[i] - numberOfCycles[mtr[r][c][id]];
                                //cerr << numberOfCycles[i] << "\t" << numberOfCycles[mtr[r][c]] << endl;
                                movesWithoutChange = 0;
                                mtr[r][c][id] = i;
                                cycleExist = (numberOfCycles[mtr[r][c][id]] > 0);
                            }
                        }
                        cycleExist = cycleExist || (numberOfCycles[mtr[r][c][id]] > 0);
                        /*if (numberOfCycles[mtr[r][c]])
                        cerr << "number of cycles = " << numberOfCycles[mtr[r][c]] << "\n\n";*/
                    }
                }

                //////////////////////////////////////////////////////////////////////////////////////////////////////
                if ((scaling) && (moves % (numberOfValuesToAssign) == 0)) {
                    numberOfActCycles = 0;

                    cycleExist = 0;
                    for (int targetId = 0; targetId < target.size(); ++targetId) {
                        ll minimum = 0;
                        for (int girth = 4; girth <= target[targetId].girth - (target[targetId].ace == 0); girth += 2) {
                            CycleEnum enumerator(girth, submatrixProto[targetId]);
                            if (!enumerator.init())
                                continue;
                            do {
                                Cycle cycle(enumerator.cycle);
                                if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                    continue;
                                processCycle(cycle, minimum, mtr, target[targetId], scalingFactor[targetId]);

                            } while (enumerator.next());
                        }

                        for (ll scale = 1; scale < maxCirc; ++scale) {
                            if (minimum == 0)
                                break;

                            ll numOfCycles = 0;
                            for (int girth = 4; girth <= target[targetId].girth - (target[targetId].ace == 0); girth += 2) {
                                CycleEnum enumerator(girth, submatrixProto[targetId]);
                                if (!enumerator.init())
                                    continue;
                                do {
                                    Cycle cycle(enumerator.cycle);
                                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                        continue;
                                    processCycle(cycle, numOfCycles, mtr, target[targetId], scale);

                                } while (enumerator.next());
                            }
                            if (numOfCycles < minimum) {
                                //numberOfActCycles += numOfCycles - minimum;

                                //cerr << numOfCycles << "\t" << minimum << endl;
                                movesWithoutChange = 0;
                                minimum = numOfCycles;
                                scalingFactor[targetId] = scale;
                            }
                        }
                        numberOfActCycles += minimum;
                        cycleExist = cycleExist || (minimum > 0);
                        //cerr << "\nsubmatrix #" << targetId + 1 << ":\nnumberOfCycles  = " << minimum << "\n\n";

                    }
                }
            }

        }
    }




    //returns number of remaining cycles
    ll anneal(bool scaling = 0, bool coupling = 0) {
        computeTotalNumberOfCycles(coupling);
        long long numberOfValuesToAssign = 0;
        for (int r = 0; r < checkNodes; ++r) {
            for (int c = 0; c < variableNodes; ++c) {
                numberOfValuesToAssign += protograph[r][c];
            }
        }
        int movesWithoutChange = 0;
        bool cycleExist = 0;
        long long moves = 0;
        //vector<ll> numCycMem;
        ll numberOfActCycles = getTotalNumberOfActiveCycles(mtr, coupling);
        ll minNumOfActiveCycles = numberOfActCycles;

        while (true) {
            for (int r = 0; r < checkNodes; ++r) {
                for (int c = 0; c < variableNodes; ++c) {
                    for (int id = 0; id < protograph[r][c]; ++id) {
                        ll curNumOfActiveCycles = numberOfActCycles;

                        ++movesWithoutChange;
                        ++moves;
                        if (numberOfActCycles == 0)
                            return 0;
                        if (stoppingCriteria(movesWithoutChange, numberOfValuesToAssign)) {
                            //noCycles = !cycleExist;
                            return curNumOfActiveCycles;

                        }
                        vector<ll> numberOfCycles(maxCirc, 0);
                        for (int targetId = 0; targetId < target.size(); ++targetId) {
                            if (r >= target[targetId].check)
                                continue;
                            if (c >= target[targetId].var)
                                continue;
                            int curGirth = target[targetId].girth - 2;
                            if (target[targetId].ace > 0)
                                curGirth += 2;
                            for (int girth = 4; girth <= curGirth; girth += 2) {
                                //cerr << r << " " << c << " " << girth << endl;
                                CycleEnum enumerator(girth, submatrixProto[targetId]);
                                if (!enumerator.init(r, c, id))
                                    continue;
                                do {
                                    Cycle cycle(enumerator.cycle);
                                    if ((enumerator.cycle[0].r != r) || (enumerator.cycle[0].c != c))
                                        break;
                                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                        continue;
                                    processCycle(cycle, numberOfCycles, mtr, target[targetId], scalingFactor[targetId]);
                                } while (enumerator.next());
                            }
                        }
                        for (int i = 0; i < maxCirc; ++i) {
                            if (i == mtr[r][c][id])
                                continue;
                            if (numberOfCycles[mtr[r][c][id]] == 0)
                                break;
                            if (numberOfCycles[i] < numberOfCycles[mtr[r][c][id]]) {
                                numberOfActCycles += numberOfCycles[i] - numberOfCycles[mtr[r][c][id]];
                                //cerr << numberOfCycles[i] << "\t" << numberOfCycles[mtr[r][c]] << endl;
                                movesWithoutChange = 0;
                                mtr[r][c][id] = i;
                                cycleExist = (numberOfCycles[mtr[r][c][id]] > 0);
                            }
                        }


                        int ind = getInd(numberOfCycles, mtr[r][c][id], moves);
                        if (numberOfCycles[ind] < numberOfCycles[mtr[r][c][id]]) {
                            movesWithoutChange = 0;
                            cycleExist = (numberOfCycles[ind] > 0);

                        }
                        numberOfActCycles += numberOfCycles[ind] - numberOfCycles[mtr[r][c][id]];
                        mtr[r][c][id] = ind;
                        cycleExist = cycleExist || (numberOfCycles[mtr[r][c][id]] > 0);
                        /*if (numberOfCycles[mtr[r][c]])
                        cerr << "number of cycles = " << numberOfCycles[mtr[r][c]] << "\n\n";*/
                    }
                }

                //////////////////////////////////////////////////////////////////////////////////////////////////////
                if ((scaling) && (moves % (5 * numberOfValuesToAssign) == 0)) {
                    numberOfActCycles = 0;

                    cycleExist = 0;
                    for (int targetId = 0; targetId < target.size(); ++targetId) {
                        ll minimum = 0;
                        for (int girth = 4; girth <= target[targetId].girth - (target[targetId].ace == 0); girth += 2) {
                            CycleEnum enumerator(girth, submatrixProto[targetId]);
                            if (!enumerator.init())
                                continue;
                            do {
                                Cycle cycle(enumerator.cycle);
                                if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                    continue;
                                processCycle(cycle, minimum, mtr, target[targetId], scalingFactor[targetId]);

                            } while (enumerator.next());
                        }

                        for (ll scale = 1; scale < maxCirc; ++scale) {
                            if (minimum == 0)
                                break;

                            ll numOfCycles = 0;
                            for (int girth = 4; girth <= target[targetId].girth - (target[targetId].ace == 0); girth += 2) {
                                CycleEnum enumerator(girth, submatrixProto[targetId]);
                                if (!enumerator.init())
                                    continue;
                                do {
                                    Cycle cycle(enumerator.cycle);
                                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                                        continue;
                                    processCycle(cycle, numOfCycles, mtr, target[targetId], scale);

                                } while (enumerator.next());
                            }
                            if (numOfCycles < minimum) {
                                //numberOfActCycles += numOfCycles - minimum;

                                //cerr << numOfCycles << "\t" << minimum << endl;
                                movesWithoutChange = 0;
                                minimum = numOfCycles;
                                scalingFactor[targetId] = scale;
                            }
                        }
                        numberOfActCycles += minimum;
                        cycleExist = cycleExist || (minimum > 0);
                        //cerr << "\nsubmatrix #" << targetId + 1 << ":\nnumberOfCycles  = " << minimum << "\n\n";
                    }
                }
            }
        }
    }


    vector<vector<vector<int> > > getMatrix() {
        return mtr;
    }
    vector<int> getScaling() {
        return scalingFactor;
    }

private:

    //if you want to change scaling procedure you need to change functions "processCycle" too
    ll getScaledValue(int x, int curCirc, int scale) {
        return ((x * scale) % maxCirc) % curCirc;

    }
    ll getTotalNumberOfActiveCycles(const vector<vector<vector<int> > >& a, bool coupling = 0) {
        ll res = 0;
        for (int targetId = 0; targetId < target.size(); ++targetId) {
            int curGirth = target[targetId].girth - 2;
            if (target[targetId].ace > 0)
                curGirth += 2;
            for (int girth = 4; girth <= curGirth; girth += 2) {
                CycleEnum enumerator(girth, submatrixProto[targetId]);
                if (!enumerator.init())
                    continue;
                do {
                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                        continue;
                    ll cur = 0;
                    bool one = true;
                    for (int i = 0; i < enumerator.cycle.size(); ++i) {
                        if (i & 1)
                            cur += getScaledValue(a[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id], target[targetId].circ, scalingFactor[targetId]);
                        else
                            cur -= getScaledValue(a[enumerator.cycle[i].r][enumerator.cycle[i].c][enumerator.cycle[i].id], target[targetId].circ, scalingFactor[targetId]);
                        if ((i) && (Tiii(enumerator.cycle[i]) == Tiii(enumerator.cycle[0])))
                            one = false;
                    }
                    if (!one) {
                        if (!lexMin(enumerator.cycle))
                            continue;
                    }
                    cur = ((cur % target[targetId].circ) + target[targetId].circ) % target[targetId].circ;
                    if (cur == 0)
                        ++res;//res += Cycle(enumerator.cycle).uniqueNodes;
                } while (enumerator.next());
            }
        }
        return res;
    }

    //get index according to distribution
    int getInd(const vector<ll>& cyc, int curInd, ll moves) {
        if (cyc[curInd] == 0)
            return curInd;
        vector<double> prob(maxCirc);
        double sumProb = 0;
        for (int i = 0; i < maxCirc; ++i) {
            prob[i] = exp(-(cyc[i] - cyc[curInd]) / temperature(totalNumberOfCycles, moves));
            /*if (i != mtr[r][c][id])
            prob[i] /= circulant;*/
            /*if ((numberOfCyclesWithTheseValues[r][c][id][i] - numberOfCyclesWithTheseValues[r][c][id][mtr[r][c][id]] == 0) && (i != mtr[r][c][id]))
            prob[i] = exp(-0.5 / temperature);*/
            sumProb += prob[i];
        }
        double randMove = sumProb * rand() / RAND_MAX;
        double sum = 0;
        for (int i = 0; i < maxCirc; ++i) {
            sum += prob[i];
            if (sum > randMove) {
                return i;
            }
        }
        return 0;
    }


    double temperature(ll totalNumberOfCycles, ll moves) {
        //cerr << 20.0 * totalNumberOfCycles / moves / moves << endl;
        return 20.0 * totalNumberOfCycles / moves / moves;
    }

    bool stoppingCriteria(int numberOfMovesWithoutChange, int numberOfValuesToAssign) {
        return numberOfMovesWithoutChange > numberOfValuesToAssign + 5;
    }
    void computeTotalNumberOfCycles(bool coupling = 0) {
        ll res = 0;
        for (int targetId = 0; targetId < target.size(); ++targetId) {
            int curGirth = target[targetId].girth - 2;
            if (target[targetId].ace > 0)
                curGirth += 2;
            for (int girth = 4; girth <= curGirth; girth += 2) {
                CycleEnum enumerator(girth, submatrixProto[targetId]);
                if (!enumerator.init())
                    continue;
                do {
                    if (!activeCycle(enumerator.cycle, target[targetId], coupling, targetId))
                        continue;
                    bool one = true;
                    for (int i = 0; i < enumerator.cycle.size(); ++i) {
                        if ((i) && (Tiii(enumerator.cycle[i]) == Tiii(enumerator.cycle[0])))
                            one = false;
                    }
                    if (!one) {
                        if (!lexMin(enumerator.cycle))
                            continue;
                    }
                    ++res;//res += Cycle(enumerator.cycle).uniqueNodes;
                } while (enumerator.next());
            }
        }
        totalNumberOfCycles = res;

    }
    bool activeCycle(const vector<entry>& cycle, const submatrix& target, bool coupling, int targetId) {
        int girth = cycle.size();
        if (girth > target.girth)
            return 0;
        int ace = 0;
        if (girth == target.girth) {
            for (int i = 0; i < cycle.size(); ++i) {
                ace += varDeg[targetId][cycle[i].c] - 2;
                //ace += checkDeg[targetId][cycle[i].r] - 2;
            }
            ace /= 2;
            if (ace >= target.ace)
                return 0;
        }
        for (int i = 0; i < cycle.size(); ++i) {
            if (cycle[i].c >= target.var)
                return 0;
            if (cycle[i].r >= target.check)
                return 0;
        }
        if (coupling) {
            long long sum = 0;
            for (int i = 0; i < cycle.size(); ++i) {
                int j = i + 1;
                if (j == cycle.size())
                    j = 0;
                if (cycle[i].r == cycle[j].r) {
                    int r = cycle[i].r;
                    int residue = target.var % target.check;
                    int shift = target.var / target.check;
                    int add = r * shift + min(r, residue);
                    int c1 = (cycle[i].c - add + target.var) % target.var;
                    int c2 = (cycle[j].c - add + target.var) % target.var;
                    sum += c2 - c1;
                }
            }
            if (sum != 0)
                return 0;
        }
        return 1;
    }
    long long gcd(long long a, long long b) {
        return b ? gcd(b, a % b) : a;
    }

    void gcd(long long a, long long b, long long& x, long long& y) {
        if (b == 0) {
            x = 1, y = 0;
            return;
        }
        gcd(b, a % b, x, y);
        long long q = a / b;
        long long xx = y;
        y = x - y * q;
        x = xx;
    }

    long long inverse(long long a, long long m) {
        long long x, y;
        gcd(a, m, x, y);
        return ((x % m) + m) % m;
    }

    bool lexMin(const Cycle& cycle) {
        int len = cycle.cycle.size();
        for (int i = 1; i < len; ++i) {
            if (cycle.cycle[i] != cycle.cycle[0])
                continue;
            if ((i & 1) == 0) {
                for (int j = i + 1, k = 1; k < len; ++j, ++k) {
                    if (j == len)
                        j = 0;
                    if (cycle.cycle[k] < cycle.cycle[j])
                        break;
                    if (cycle.cycle[j] < cycle.cycle[k])
                        return false;

                }
            }
            else {
                for (int j = i - 1, k = 1; k < len; --j, ++k) {
                    if (j == -1)
                        j = len - 1;
                    if (cycle.cycle[j] < cycle.cycle[k])
                        return false;
                    if (cycle.cycle[k] < cycle.cycle[j])
                        break;
                }
            }
        }
        return true;
    }

    void processCycle(const Cycle& cycle, vector<ll>& numberOfCycles, const vector<vector<vector<int> > >& a, const submatrix& target, int scalingFactor, ll mult = 1) {
        mult = cycle.uniqueNodes;
        Tiii cur = cycle.cycle[0];
        long long c = 0, b = 0;//cx=b mod circulant
        bool one = true;
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            if (cur == cycle.cycle[i]) {
                if (i & 1)
                    ++c;
                else
                    --c;
                if (i)
                    one = false;
            }
            else {
                if (i & 1)
                    b -= getScaledValue(a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third], target.circ, scalingFactor);
                else
                    b += getScaledValue(a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third], target.circ, scalingFactor);

            }
        }
        if (!one) {
            if (!lexMin(cycle))
                return;
        }
        c = ((c % target.circ) + target.circ) % target.circ;
        b = ((b % target.circ) + target.circ) % target.circ;
        if (c == 0) {
            if (b != 0)
                return;
            for (int i = 0; i < maxCirc; ++i)
                numberOfCycles[i] += mult;
            return;
        }
        long long d1 = gcd(c, target.circ);
        if ((b % d1) != 0)
            return;
        int circ = target.circ / d1;
        c /= d1, b /= d1;
        int invr = inverse(scalingFactor, maxCirc);
        long long x = (inverse(c, circ) * b) % circ;
        for (int i = 0; x + circ * i < target.circ; ++i) {
            int u = x + circ * i;
            for (int j = 0; u + j * target.circ < maxCirc; ++j) {
                numberOfCycles[((u + j * target.circ) * invr) % maxCirc] += mult;

            }
        }
    }


    void processCycle(const Cycle& cycle, ll& numberOfCycles, const vector<vector<vector<int> > >& a, const submatrix& target, ll scalingFactor, ll mult = 1) {
        //mult = cycle.uniqueNodes;
        Tiii cur = cycle.cycle[0];
        bool one = true;
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            if (cur == cycle.cycle[i]) {
                if (i)
                    one = false;
            }
        }
        if (!one) {
            if (!lexMin(cycle))
                return;
        }
        long long b = 0;
        //
        for (int i = 0; i < cycle.cycle.size(); ++i) {
            if (i & 1)
                b -= getScaledValue(a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third], target.circ, scalingFactor);
            else
                b += getScaledValue(a[cycle.cycle[i].first][cycle.cycle[i].second][cycle.cycle[i].third], target.circ, scalingFactor);
        }
        b = ((b % target.circ) + target.circ) % target.circ;
        if (b == 0)
            numberOfCycles += mult;

    }
};






