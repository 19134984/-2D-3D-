#include <unordered_map>

template <typename T>
class ReferenceData{

    private:
    T compVals[4][6] = {
      { 164.24, 701.92, 4.831, 0.851, 0.020, 16.790       },
      { 389.88, 2241.37, 5.749, 0.937, 0.011, 30.506      },
      { 503.24, 6820.07, 13.552, 0.966, 0.0064, 57.350    },
      { 2323.00, 21463.00, 9.239, 0.940, 0.491, 103.663   }
    };

    public:
    ReferenceData(){};

    bool hasRayleigh(std::size_t Ra) {
        std::unordered_map<T, T> mult_map =
            {
                {1e7,   true},
                {1e8,   true},
                {1e9,   true},
                {1e10,  true}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : false;
    }

    T multiplyCharU(std::size_t Ra, T charU){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,   compVals[0][1]},
                {1e8,   compVals[1][1]},
                {1e9,   compVals[2][1]},
                {1e10,  compVals[3][1]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] * charU : charU;
    }

    T getUx_max(std::size_t Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,  compVals[0][0]},
                {1e8,  compVals[1][0]},
                {1e9,  compVals[2][0]},
                {1e10,  compVals[3][0]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getUy_max(std::size_t Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,  compVals[0][1]},
                {1e8,  compVals[1][1]},
                {1e9,  compVals[2][1]},
                {1e10,  compVals[3][1]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getY_max(std::size_t Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,  compVals[0][3]},
                {1e8,  compVals[1][3]},
                {1e9,  compVals[2][3]},
                {1e10,  compVals[3][3]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getX_max(std::size_t Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,  compVals[0][4]},
                {1e8,  compVals[1][4]},
                {1e9,  compVals[2][4]},
                {1e10,  compVals[3][4]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

    T getNusselt(std::size_t Ra){
        std::unordered_map<T, T> mult_map =
            {
                {1e7,  compVals[0][5]},
                {1e8,  compVals[1][5]},
                {1e9,  compVals[2][5]},
                {1e10,  compVals[3][5]}
            };

        return mult_map.contains(Ra) ? mult_map[Ra] : -1;
    }

};
