/*
 * ============================================================================
 * MÉTODO EXPLÍCITO PARA LA ECUACIÓN DEL CALOR 1D
 * Ecuación: ∂u/∂t = α² ∂²u/∂x²
 * Esquema: u_i^(n+1) = λ·u_(i-1)^n + (1-2λ)·u_i^n + λ·u_(i+1)^n
 * Condición de estabilidad: λ ≤ 0.5
 * ============================================================================
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <string>
#include <algorithm>

using namespace std;

const double PI = 3.141592653589793;

// ============================================================================
// CLASE: HeatEquationExplicit
// ============================================================================

class HeatEquationExplicit {
private:
    // Parámetros físicos
    double L;       // Longitud del dominio
    double alpha;   // Coeficiente de difusión
    double T;       // Tiempo final
    
    // Parámetros numéricos
    int m;          // Divisiones espaciales
    int nt;         // Pasos temporales
    double h;       // Paso espacial
    double k;       // Paso temporal
    double lambda;  // Parámetro λ = α²k/h²
    
    // Arrays
    vector<double> x;      // Nodos espaciales
    vector<double> u;      // Solución actual
    vector<double> u_new;  // Solución en el siguiente paso
    
    // Métricas
    double execution_time;
    double error_max;
    double error_l2;
    
public:
    // Constructor
    HeatEquationExplicit(double L_val = 1.0, double alpha_val = 1.0, 
                         double T_val = 0.1, int m_val = 40, double lambda_val = 0.4)
        : L(L_val), alpha(alpha_val), T(T_val), m(m_val), lambda(lambda_val) {
        
        // Calcular parámetros
        h = L / m;
        k = lambda * (h * h) / (alpha * alpha);
        nt = static_cast<int>(ceil(T / k));
        k = T / nt;  // Ajustar para cubrir exactamente T
        lambda = alpha * alpha * k / (h * h);  // Recalcular lambda
        
        // Inicializar arrays
        x.resize(m + 1);
        u.resize(m + 1);
        u_new.resize(m + 1);
        
        // Condición inicial
        initialize();
    }
    
    // Inicializar condición inicial y malla
    void initialize() {
        for (int i = 0; i <= m; ++i) {
            x[i] = i * h;
            u[i] = sin(PI * x[i] / L);
        }
        // Condiciones de frontera
        u[0] = 0.0;
        u[m] = 0.0;
    }
    
    // Solución exacta
    double exact_solution(double x_val, double t) const {
        return sin(PI * x_val / L) * exp(-alpha * alpha * PI * PI * t / (L * L));
    }
    
    // Un paso del método explícito
    void explicit_step() {
        // Puntos interiores
        for (int i = 1; i < m; ++i) {
            u_new[i] = lambda * u[i-1] + (1.0 - 2.0 * lambda) * u[i] + lambda * u[i+1];
        }
        // Condiciones de frontera
        u_new[0] = 0.0;
        u_new[m] = 0.0;
        
        // Actualizar solución
        u = u_new;
    }
    
    // Resolver la ecuación completa
    void solve() {
        auto start = chrono::high_resolution_clock::now();
        
        for (int n = 0; n < nt; ++n) {
            explicit_step();
        }
        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = end - start;
        execution_time = diff.count();
    }
    
    // Calcular errores
    void compute_errors() {
        error_max = 0.0;
        error_l2 = 0.0;
        
        for (int i = 0; i <= m; ++i) {
            double u_exact = exact_solution(x[i], T);
            double diff = fabs(u[i] - u_exact);
            error_max = max(error_max, diff);
            error_l2 += diff * diff;
        }
        
        error_l2 = sqrt(error_l2 * h);
    }
    
    // Verificar estabilidad
    bool is_stable() const {
        return lambda <= 0.5;
    }
    
    // Getters
    int get_m() const { return m; }
    int get_nt() const { return nt; }
    double get_h() const { return h; }
    double get_k() const { return k; }
    double get_lambda() const { return lambda; }
    double get_time() const { return execution_time; }
    double get_error_max() const { return error_max; }
    double get_error_l2() const { return error_l2; }
    const vector<double>& get_solution() const { return u; }
    const vector<double>& get_x() const { return x; }
};

// ============================================================================
// FUNCIONES DE EXPERIMENTOS
// ============================================================================

void print_header(const string& title) {
    cout << "\n" << string(80, '=') << endl;
    cout << title << endl;
    cout << string(80, '-') << endl;
}

void experiment_stability() {
    print_header("EXPERIMENTO 1: ANÁLISIS DE ESTABILIDAD");
    cout << "Condición de estabilidad: λ ≤ 0.5\n" << endl;
    
    vector<double> lambda_values = {0.3, 0.45, 0.5, 0.55, 0.75, 1.0};
    int m = 40;
    
    cout << left << setw(10) << "Lambda" 
         << setw(12) << "Estable?" 
         << setw(12) << "h"
         << setw(14) << "k"
         << setw(10) << "nt"
         << setw(15) << "Error Max"
         << setw(13) << "Error L2"
         << setw(12) << "Tiempo(s)" << endl;
    cout << string(98, '-') << endl;
    
    for (double lam : lambda_values) {
        HeatEquationExplicit solver(1.0, 1.0, 0.1, m, lam);
        solver.solve();
        solver.compute_errors();
        
        cout << fixed << setprecision(3) << setw(10) << lam
             << setw(12) << (solver.is_stable() ? "SI" : "NO")
             << setprecision(6) << setw(12) << solver.get_h()
             << setprecision(8) << setw(14) << solver.get_k()
             << setw(10) << solver.get_nt()
             << scientific << setprecision(6) << setw(15) << solver.get_error_max()
             << setw(13) << solver.get_error_l2()
             << fixed << setprecision(6) << setw(12) << solver.get_time() << endl;
    }
    
    cout << "\nOBSERVACIONES:" << endl;
    cout << "  - Para λ ≤ 0.5: El método es estable y convergente" << endl;
    cout << "  - Para λ > 0.5: El método es INESTABLE (errores crecen)" << endl;
}

void experiment_convergence() {
    print_header("EXPERIMENTO 2: ESTUDIO DE CONVERGENCIA");
    cout << "λ fijo = 0.4 (estable)\n" << endl;
    
    vector<int> m_values = {10, 20, 40, 80, 160, 320};
    double lambda = 0.4;
    
    cout << left << setw(8) << "m"
         << setw(12) << "h"
         << setw(14) << "k"
         << setw(10) << "nt"
         << setw(15) << "Error Max"
         << setw(13) << "Error L2"
         << setw(12) << "Tiempo(s)"
         << setw(15) << "Nodos·Pasos" << endl;
    cout << string(99, '-') << endl;
    
    for (int m : m_values) {
        HeatEquationExplicit solver(1.0, 1.0, 0.1, m, lambda);
        solver.solve();
        solver.compute_errors();
        
        cout << setw(8) << m
             << fixed << setprecision(6) << setw(12) << solver.get_h()
             << setprecision(8) << setw(14) << solver.get_k()
             << setw(10) << solver.get_nt()
             << scientific << setprecision(6) << setw(15) << solver.get_error_max()
             << setw(13) << solver.get_error_l2()
             << fixed << setprecision(6) << setw(12) << solver.get_time()
             << setw(15) << (m+1) * solver.get_nt() << endl;
    }
    
    cout << "\nOBSERVACIONES:" << endl;
    cout << "  - Al reducir h, el error disminuye como O(h²)" << endl;
    cout << "  - El número de pasos temporales crece como nt ~ 1/h²" << endl;
    cout << "  - Costo computacional total ~ O(1/h³)" << endl;
}

void experiment_performance() {
    print_header("EXPERIMENTO 3: ANÁLISIS DE RENDIMIENTO");
    cout << "Probando mallas muy refinadas (λ = 0.45)\n" << endl;
    
    vector<int> m_values = {100, 200, 500, 1000, 1500, 2000, 3000, 4000};
    double lambda = 0.45;
    
    cout << left << setw(10) << "m"
         << setw(12) << "Nodos"
         << setw(10) << "nt"
         << setw(14) << "k"
         << setw(13) << "Ops·10⁶"
         << setw(12) << "Tiempo(s)"
         << setw(12) << "Mem(KB)" << endl;
    cout << string(83, '-') << endl;
    
    for (int m : m_values) {
        HeatEquationExplicit solver(1.0, 1.0, 0.1, m, lambda);
        solver.solve();
        
        long long operations = (long long)(m - 1) * solver.get_nt() * 5;
        double ops_millions = operations / 1.0e6;
        double memory_kb = (m + 1) * sizeof(double) * 3 / 1024.0;
        
        cout << setw(10) << m
             << setw(12) << (m + 1)
             << setw(10) << solver.get_nt()
             << fixed << setprecision(8) << setw(14) << solver.get_k()
             << setprecision(3) << setw(13) << ops_millions
             << setprecision(6) << setw(12) << solver.get_time()
             << setprecision(2) << setw(12) << memory_kb << endl;
    }
    
    cout << "\nLÍMITE PRÁCTICO:" << endl;
    cout << "  - Para m ~ 2000: nt ~ 1.6·10⁶ pasos temporales" << endl;
    cout << "  - Costo: ~ 16 millones de operaciones" << endl;
    cout << "  - Limitación principal: Restricción de estabilidad k ~ h²" << endl;
}

void experiment_comparison() {
    print_header("EXPERIMENTO 4: COMPARACIÓN DE MÉTODOS");
    
    cout << "\nCaracterísticas del método EXPLÍCITO:" << endl;
    cout << "  ✓ Implementación muy simple" << endl;
    cout << "  ✓ No requiere resolver sistemas de ecuaciones" << endl;
    cout << "  ✓ Bajo costo por paso temporal: O(m)" << endl;
    cout << "  ✗ CONDICIONALMENTE ESTABLE (λ ≤ 0.5)" << endl;
    cout << "  ✗ Requiere muchos pasos temporales para mallas finas" << endl;
    cout << "  ✗ Error temporal O(k)" << endl;
    
    cout << "\nCaracterísticas del método CRANK-NICOLSON (Implícito):" << endl;
    cout << "  ✓ Incondicionalmente estable (cualquier λ)" << endl;
    cout << "  ✓ Mayor precisión temporal O(k²)" << endl;
    cout << "  ✓ Permite pasos temporales más grandes" << endl;
    cout << "  ✗ Requiere resolver sistemas tridiagonales: O(m)" << endl;
    cout << "  ✗ Implementación más compleja" << endl;
    
    cout << "\nCOMPARACIÓN DIRECTA (m = 80):" << endl;
    cout << string(80, '-') << endl;
    
    int m = 80;
    
    // Método explícito con λ = 0.45
    HeatEquationExplicit explicit_solver(1.0, 1.0, 0.1, m, 0.45);
    explicit_solver.solve();
    explicit_solver.compute_errors();
    
    cout << "\nMétodo EXPLÍCITO (λ = 0.45):" << endl;
    cout << "  Pasos temporales: " << explicit_solver.get_nt() << endl;
    cout << "  k = " << scientific << setprecision(6) << explicit_solver.get_k() << endl;
    cout << "  Error máximo: " << explicit_solver.get_error_max() << endl;
    cout << "  Error L2: " << explicit_solver.get_error_l2() << endl;
    cout << "  Tiempo: " << fixed << setprecision(6) << explicit_solver.get_time() << " s" << endl;
    
    cout << "\nMétodo CRANK-NICOLSON (λ = 2.0, tu implementación actual):" << endl;
    cout << "  Podría usar k ~ 4.4 veces mayor" << endl;
    cout << "  Requeriría ~4 veces menos pasos temporales" << endl;
    cout << "  Pero cada paso cuesta más (resolver sistema tridiagonal)" << endl;
}

// ============================================================================
// PROGRAMA PRINCIPAL
// ============================================================================

int main() {
    cout << string(80, '=') << endl;
    cout << "  MÉTODO DE DIFERENCIAS FINITAS EXPLÍCITO - ECUACIÓN DEL CALOR 1D" << endl;
    cout << string(80, '=') << endl;
    
    // Ejemplo básico
    print_header("EJEMPLO BÁSICO");
    
    HeatEquationExplicit solver(1.0, 1.0, 0.1, 40, 0.4);
    
    cout << "\nParámetros:" << endl;
    cout << "  L = 1.0, α = 1.0, T = 0.1" << endl;
    cout << "  m = " << solver.get_m() << ", λ = " << fixed << setprecision(4) 
         << solver.get_lambda() << endl;
    cout << "  h = " << setprecision(6) << solver.get_h() 
         << ", k = " << setprecision(8) << solver.get_k() << endl;
    cout << "  Estabilidad: " << (solver.is_stable() ? "✓ ESTABLE" : "✗ INESTABLE") << endl;
    
    solver.solve();
    solver.compute_errors();
    
    cout << "\nResultados:" << endl;
    cout << "  Error máximo: " << scientific << setprecision(6) 
         << solver.get_error_max() << endl;
    cout << "  Error L2: " << solver.get_error_l2() << endl;
    cout << "  Tiempo de ejecución: " << fixed << setprecision(6) 
         << solver.get_time() << " s" << endl;
    
    // Mostrar algunos valores de la solución
    cout << "\nPrimeros valores de la solución en t = T:" << endl;
    const auto& u = solver.get_solution();
    const auto& x = solver.get_x();
    for (int i = 0; i <= min(10, solver.get_m()); ++i) {
        double u_exact = solver.exact_solution(x[i], 0.1);
        cout << "  x[" << setw(2) << i << "] = " << fixed << setprecision(4) 
             << x[i] << ", u = " << setprecision(6) << u[i] 
             << ", u_exact = " << u_exact 
             << ", error = " << scientific << setprecision(2) 
             << fabs(u[i] - u_exact) << endl;
    }
    
    // Experimentos
    experiment_stability();
    experiment_convergence();
    experiment_performance();
    experiment_comparison();
    
    // Resumen final
    print_header("RESUMEN DE CONCLUSIONES");
    
    cout << "\n1. ESTABILIDAD:" << endl;
    cout << "   - El método explícito es estable SOLO si λ ≤ 0.5" << endl;
    cout << "   - Para λ > 0.5, las soluciones divergen (inestabilidad numérica)" << endl;
    cout << "   - Condición: k ≤ h²/(2α²)" << endl;
    
    cout << "\n2. CONVERGENCIA:" << endl;
    cout << "   - Error espacial: O(h²)" << endl;
    cout << "   - Error temporal: O(k)" << endl;
    cout << "   - Error global: O(h² + k)" << endl;
    
    cout << "\n3. COSTO COMPUTACIONAL:" << endl;
    cout << "   - Operaciones por paso: 5(m-1) ~ O(m)" << endl;
    cout << "   - Pasos temporales: nt ~ T/(λh²/α²) ~ O(1/h²)" << endl;
    cout << "   - Costo total: O(m·nt) = O(m/h²) = O(1/h³)" << endl;
    
    cout << "\n4. COMPARACIÓN CON CRANK-NICOLSON:" << endl;
    cout << "   - Explícito: Simple pero limitado por estabilidad" << endl;
    cout << "   - Crank-Nicolson: Más complejo pero sin restricción en k" << endl;
    cout << "   - Para mallas finas: Crank-Nicolson es más eficiente" << endl;
    
    cout << "\n5. RECOMENDACIONES:" << endl;
    cout << "   - Usar método explícito para: prototipos, mallas gruesas" << endl;
    cout << "   - Usar método implícito para: producción, mallas finas" << endl;
    cout << "   - Siempre verificar λ ≤ 0.5 en métodos explícitos" << endl;
    
    cout << "\n" << string(80, '=') << endl;
    cout << "ANÁLISIS COMPLETADO" << endl;
    cout << string(80, '=') << endl << endl;
    
    return 0;
}
