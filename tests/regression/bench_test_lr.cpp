#include <iostream>
#include <chrono>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace boost::math;
using boost::multiprecision::cpp_dec_float_50;

template <typename T>
double compute_p_value(T F_stat, T df_num, T df_den)
{
    fisher_f_distribution<T> dist(df_num, df_den);
    T p = T(1) - cdf(dist, F_stat);
    return static_cast<double>(p);
}

int main()
{
    const int N = 1000;
    double df_num = 5.0, df_den = 200.0;

    // Dataset 1: non-significant F-stat
    double F_ns = 1.2;
    // Dataset 2: highly significant F-stat
    double F_sig = 120.0;

    // --- Double precision ---
    {
        auto start = chrono::high_resolution_clock::now();
        double p = 0.0;
        for (int i = 0; i < N; ++i)
        {
            p = compute_p_value(F_ns, df_num, df_den);
            p = compute_p_value(F_sig, df_num, df_den);
        }
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "[double] time: " << elapsed.count() << " s\n";
        cout << "  p(non-sig) = " << compute_p_value(F_ns, df_num, df_den) << "\n";
        cout << "  p(sig)     = " << compute_p_value(F_sig, df_num, df_den) << "\n";
    }

    // --- High precision (cpp_dec_float_50) ---
    {
        cpp_dec_float_50 F_ns_hp = F_ns;
        cpp_dec_float_50 F_sig_hp = F_sig;
        cpp_dec_float_50 df_num_hp = df_num;
        cpp_dec_float_50 df_den_hp = df_den;

        auto start = chrono::high_resolution_clock::now();
        double p = 0.0;
        for (int i = 0; i < N; ++i)
        {
            p = compute_p_value(F_ns_hp, df_num_hp, df_den_hp);
            p = compute_p_value(F_sig_hp, df_num_hp, df_den_hp);
        }
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = end - start;
        cout << "[cpp_dec_float_50] time: " << elapsed.count() << " s\n";
        cout << "  p(non-sig) = " << compute_p_value(F_ns_hp, df_num_hp, df_den_hp) << "\n";
        cout << "  p(sig)     = " << compute_p_value(F_sig_hp, df_num_hp, df_den_hp) << "\n";
    }

    return 0;
}
