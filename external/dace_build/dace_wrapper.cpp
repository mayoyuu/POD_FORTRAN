// 文件名: dace_wrapper.cpp
#include <dace/dace.h>
#include <vector>
#include <iostream>
#include <cstring>

using namespace DACE;

// 智能对象池与回收站
static std::vector<DA> da_registry;
static std::vector<int> free_handles;
// 在文件顶部全局变量区补充 compiledDA 的智能对象池
static std::vector<compiledDA> cda_registry;
static std::vector<int> free_cda_handles;

extern "C" {
    // ==========================================
    // DACE 引擎全局控制 (必须走 C++ 的静态方法)
    // ==========================================
    void fdace_core_init(int ord, int nvar) {
        // 调用 C++ 类的静态初始化方法，它会自动设置 initialized 标志并调用底层 daceInitialize
        DACE::DA::init(ord, nvar); 
    }

    // ==========================================
    // 内存与生命周期管理
    // ==========================================
    void fdace_allocate(int* handle) {
        if (!free_handles.empty()) {
            *handle = free_handles.back();
            free_handles.pop_back();
            da_registry[*handle] = DA(0.0); // 赋初值0
        } else {
            da_registry.push_back(DA(0.0));
            *handle = da_registry.size() - 1;
        }
    }

    // 初始化为独立变量 (如 x1, x2)
    void fdace_allocate_var(int* handle, int var_idx) {
        fdace_allocate(handle);
        da_registry[*handle] = DA(var_idx, 1.0); 
    }

    // 释放句柄 (极其重要！防止内存泄漏)
    void fdace_free(int handle) {
        if (handle >= 0 && handle < da_registry.size()) {
            da_registry[handle] = DA(0.0); // 释放 DACE 底层内存
            free_handles.push_back(handle);
        }
    }

    void fdace_copy(int h_src, int h_dest) {
        da_registry[h_dest] = da_registry[h_src];
    }
    // ==========================================
    // compiledDA 内存与生命周期管理
    // ==========================================
    
    // 将整个 DA 数组（例如状态向量）一次性编译为一个 compiledDA 对象
    void fdace_compile_vector(const int* da_handles, int size, int* cda_handle) {
        std::vector<DA> vec(size);
        for(int i = 0; i < size; ++i) {
            vec[i] = da_registry[da_handles[i]];
        }
        
        // 核心：调用 DACE 的批量编译构造函数
        compiledDA cda(vec); 
        
        if (!free_cda_handles.empty()) {
            *cda_handle = free_cda_handles.back();
            free_cda_handles.pop_back();
            cda_registry[*cda_handle] = cda;
        } else {
            cda_registry.push_back(cda);
            *cda_handle = cda_registry.size() - 1;
        }
    }

    // 释放 compiledDA 句柄
    void fdace_compiled_free(int handle) {
        if (handle >= 0 && handle < cda_registry.size()) {
            // 用一个空的 DA 的编译结果覆盖，释放底层 ac 数组内存
            cda_registry[handle] = compiledDA(DA(0.0)); 
            free_cda_handles.push_back(handle);
        }
    }

    // ==========================================
    // 极速模板求值接口 (显式转换为 std::vector<double>)
    // ==========================================
    void fdace_compiled_eval_double(int cda_handle, const double* in_args, int n_args, double* out_res, int n_res) {
        // 利用 std::vector 包装输入，以触发 DACE 内部最高效的 double 特化求值
        std::vector<double> args_vec(in_args, in_args + n_args);
        std::vector<double> res_vec(n_res);
        
        // 核心求值：这一步的速度极快
        cda_registry[cda_handle].eval(args_vec, res_vec); 
        
        // 将结果写回 C 数组，供 Fortran 读取
        for(int i = 0; i < n_res; ++i) {
            out_res[i] = res_vec[i];
        }
    }

    // ==========================================
    // 基础代数运算
    // ==========================================
    // 实数赋值 (覆盖原有内存中的多项式)
    void fdace_add(int h1, int h2, int h_out) { da_registry[h_out] = da_registry[h1] + da_registry[h2]; }
    void fdace_sub(int h1, int h2, int h_out) { da_registry[h_out] = da_registry[h1] - da_registry[h2]; }
    void fdace_mul(int h1, int h2, int h_out) { da_registry[h_out] = da_registry[h1] * da_registry[h2]; }
    void fdace_div(int h1, int h2, int h_out) { da_registry[h_out] = da_registry[h1] / da_registry[h2]; }

    void fdace_add_double(int h1, double val, int h_out) { da_registry[h_out] = da_registry[h1] + val; }
    void fdace_mul_double(int h1, double val, int h_out) { da_registry[h_out] = da_registry[h1] * val; }

    // ==========================================
    // 补充的实数与一元算符接口
    // ==========================================
    // 一元负号 (-DA)
    void fdace_negate(int h_in, int h_out) {da_registry[h_out] = -da_registry[h_in];}
    // DA - 实数
    void fdace_sub_double(int h1, double val, int ho) {da_registry[ho] = da_registry[h1] - val;}
    // 实数 - DA
    void fdace_double_sub(double val, int h2, int ho) {da_registry[ho] = val - da_registry[h2];}
    // 实数 / DA
    void fdace_double_div(double val, int h2, int ho) {da_registry[ho] = val / da_registry[h2];}

    // ==========================================
    // 幂次运算 (Power)
    // ==========================================
    
    // DA 的整数次幂: out = in ^ p
    void fdace_pow_int(int h_in, int p, int h_out) {
        // DACE 重载了标准库的 pow
        da_registry[h_out] = pow(da_registry[h_in], p);
    }

    // DA 的实数次幂: out = in ^ p
    void fdace_pow_double(int h_in, double p, int h_out) {
        da_registry[h_out] = pow(da_registry[h_in], p);
    }

    // ==========================================
    // 科学与微积分函数
    // ==========================================
    void fdace_sin(int h_in, int h_out) { da_registry[h_out] = sin(da_registry[h_in]); }
    void fdace_cos(int h_in, int h_out) { da_registry[h_out] = cos(da_registry[h_in]); }
    void fdace_exp(int h_in, int h_out) { da_registry[h_out] = exp(da_registry[h_in]); }
    void fdace_sqrt(int h_in, int h_out) { da_registry[h_out] = sqrt(da_registry[h_in]); }

    void fdace_asin(int h_in, int h_out) { 
        da_registry[h_out] = asin(da_registry[h_in]); 
    }

    void fdace_atan2(int hy, int hx, int h_out) { 
        da_registry[h_out] = atan2(da_registry[hy], da_registry[hx]); 
    }
    
    // 求导 (对第 var_idx 个变量求导)
    void fdace_deriv(int h_in, int var_idx, int h_out) { 
        da_registry[h_out] = da_registry[h_in].deriv(var_idx); 
    }

   // 系数提取器 (纯净透传)
    double c_fdace_get_coeff(int handle, const unsigned int* exponents) {
        const unsigned int nvar = daceGetMaxVariables();
        std::vector<unsigned int> jj(exponents, exponents + nvar);
        return da_registry[handle].getCoefficient(jj);
    }

    // 获取当前引擎的最大变量数
    int c_fdace_get_max_variables() {
        return daceGetMaxVariables();
    }

    // 获取常数项
    void fdace_get_cons(int handle, double* val) {
        *val = da_registry[handle].cons();
    }

    void fdace_print(int handle) {
        std::cout << da_registry[handle] << std::endl;
    }

    // ==========================================
    // 代入与求值 (Evaluation / Plug)
    // ==========================================
    // 将 DA 多项式中的第 var_idx 个变量替换为实数值 val
    void fdace_eval_var(int h_in, int var_idx, double val, int h_out) {
        da_registry[h_out] = da_registry[h_in].plug(var_idx, val);
    }

    // ==========================================
    // 全代入求值 (Full Evaluation)
    // ==========================================
    // 传入一个实数数组，将 DA 中所有的变量一次性全部替换，返回最终的常数值
    // void fdace_eval_all(int h_in, const double* vals, int n_vals, double* result) {
    //     DA temp = da_registry[h_in];
    //     // 循环调用 plug，将所有变量依次代入
    //     for (int i = 0; i < n_vals; ++i) {
    //         temp = temp.plug(i + 1, vals[i]); // DACE 变量索引从 1 开始
    //     }
    //     // 提取最终的常数项并返回
    //     *result = temp.cons();
    // }
    // ==========================================
    // 全代入求值 (极速版 - 替换原有基于 plug 循环的版本)
    // ==========================================
    void fdace_eval_all(int h_in, const double* vals, int n_vals, double* result) {
        // 将传入的 C 数组包装为 std::vector<double>
        std::vector<double> args_vec(vals, vals + n_vals);
        
        // 直接调用 DACE 官方的 eval 模板。
        // 这一步在 DACE 内部会自动：
        // 1. 生成临时的 compiledDA
        // 2. 调用高度优化的 double 批量求值机制
        // 3. 返回一个单独的 double
        // 4. 自动销毁临时的 compiledDA
        *result = da_registry[h_in].eval(args_vec);
    }

    // ==========================================
    // 截断与阶数控制 (Truncation & Order Control)
    // ==========================================
    
    // 1. 全局截断阶数控制
    void fdace_set_to(int ot) { DACE::DA::setTO(ot); }
    void fdace_push_to(int ot) { DACE::DA::pushTO(ot); }
    void fdace_pop_to() { DACE::DA::popTO(); }

    // 2. 基于数值大小的微小系数截断 (Epsilon)
    void fdace_set_eps(double eps) { DACE::DA::setEps(eps); }

    // 3. 特定多项式的阶数裁剪 (Trimming)
    void fdace_trim(int h_in, int min_ord, int max_ord, int h_out) {
        da_registry[h_out] = da_registry[h_in].trim(min_ord, max_ord);
    }

    // 4. 常数项的数学取整截断 (Trunc)
    void fdace_trunc(int h_in, int h_out) {
        da_registry[h_out] = da_registry[h_in].trunc();
    }

    // 获取当前截断阶数
    int fdace_get_to() { return DACE::DA::getTO(); }
}