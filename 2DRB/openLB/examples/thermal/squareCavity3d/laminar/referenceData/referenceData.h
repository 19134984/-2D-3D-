#include <unordered_map>

template <typename T>
class ReferenceData{

    private:
    T compVals[4][6] = {
            {  3.649,  3.696, 1.013, 0.813, 0.178, 1.117 },
            { 16.178, 19.617, 1.212, 0.823, 0.119, 2.238 },
            { 34.730, 68.590, 1.975, 0.855, 0.066, 4.509 },
            { 64.530, 219.36, 3.400, 0.850, 0.036, 8.817 }
        };

    public:
    ReferenceData(){};

    bool hasRayleigh(T Ra) {
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  true},
                {1e4,  true},
                {1e5,  true},
                {1e6,  true}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : false;
    }

    T multiplyCharU(T Ra, T charU){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][1]},
                {1e4,  compVals[1][1]},
                {1e5,  compVals[2][1]},
                {1e6,  compVals[3][1]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] * charU : charU;
    }

    T getUx_max(T Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][0]},
                {1e4,  compVals[1][0]},
                {1e5,  compVals[2][0]},
                {1e6,  compVals[3][0]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getUy_max(T Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][1]},
                {1e4,  compVals[1][1]},
                {1e5,  compVals[2][1]},
                {1e6,  compVals[3][1]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getY_max(T Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][3]},
                {1e4,  compVals[1][3]},
                {1e5,  compVals[2][3]},
                {1e6,  compVals[3][3]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getX_max(T Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][4]},
                {1e4,  compVals[1][4]},
                {1e5,  compVals[2][4]},
                {1e6,  compVals[3][4]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getNusselt(T Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e3,  compVals[0][5]},
                {1e4,  compVals[1][5]},
                {1e5,  compVals[2][5]},
                {1e6,  compVals[3][5]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

};
