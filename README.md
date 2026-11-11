# hamiltonian_systems

**hamiltonian_systems** — header-only тулкит на **C++23** для интегрирования гамильтоновых систем; основной акцент — **небесная механика**. Проект объединяет расщеплённые (симплектические) методы —  **Leapfrog** и семейство Йошиды — со схемами Рунге-Кутты, а также **символьным дифференцированием в момент компиляции**: из аналитического выражения гамильтониана автоматически строится правая часть канонических уравнений.

## Содержание

- [Ключевые возможности](#ключевые-возможности)
- [Структура репозитория](#структура-репозитория)
- [Требования и сборка](#требования-и-сборка)
- [Быстрый старт](#быстрый-старт)
- **Модули**
  - [precision](#precision)
  - [auto_diff](#auto_diff)
  - [Гамильтоновы системы (ham.hpp)](#гамильтоновы-системы-hamhpp)
  - [newton](#newton)
  - [runge_kutta](#runge_kutta)
  - [parallel](#parallel)
  - [yoshida](#yoshida)
- [Примеры](#примеры)
- [Тесты](#тесты)
- [Лицензия и стиль](#лицензия-и-стиль)

---

## Ключевые возможности

* **Симметричные симплектические методы** — `Leapfrog` и семейство **Yoshida** (чётные порядки $2, 4, 6, 8, \dots$), построенные для расщеплённых гамильтонианов.
* **Символьное дифференцирование в момент компиляции** (`auto_diff`): выражение представляется *типом*, и производные вычисляются на этапе компиляции.
* **Автоматическая сборка правой части** (`ham.hpp`): `Hamiltonian::compile_system<N>` превращает выражение гамильтониана в RHS-функтор канонических уравнений.
* **Методы Рунге-Кутты** — классический `RK4` и симплектический коллокационный `BRK` (узлы сетки Гаусса-Лежандра),.
* **Параллелизация OpenMP** — `for_i`, `for_each`, `reduce`.

---

## Структура репозитория

```
hamiltonian_systems/
├── precision/        # тип Real и обёртки mp_*
├── auto_diff/        # символьное дифференцирование (Expression-дерева)
│   ├── auto_diff.hpp
│   ├── auto_diff_base.hpp
│   ├── auto_diff_expressions.hpp
│   ├── auto_diff_normalizer.hpp
│   ├── auto_diff_squise.hpp
│   ├── auto_diff_subst.hpp
│   ├── auto_diff_squiser_rules.hpp
│   └── ham.hpp       # Hamiltonian::compile_system
├── newton/           # метод Ньютона, линейные решатели, минимизация
│   ├── newton.hpp
│   └── src/newton.cpp
├── runge_kutta/      # RK4, BRK, Integrate
│   └── include/runge_kutta/{odesolver,rk4,brk}.hpp
├── parallel/         # OpenMP-хелперы
│   └── parallel.hpp
├── yoshida/          # Yoshida и Leapfrog
│   └── include/yoshida/{yoshida,leapfrog,leapfrog_extension}.hpp
├── examples/         # pendulum, sitnikov, three_body
├── tests/            # doctest-тесты (auto_diff)
└── CMakeLists.txt
```

---

## Требования и сборка

**Зависимости:**
* CMake $\ge 3.25$;
* компилятор с поддержкой **C++23** (по умолчанию `clang`/`clang++`);
* **Boost** (в т.ч. `boost::numeric::ublas`, `boost::multiprecision`);
* **OpenMP**;
* **GMP/MPFR** — только при `USE_MULTIPRECISION`;
* **doctest** — только для тестов.

**Сборка:**
```bash
cmake -S . -B build -G Ninja
cmake --build build
```

**Высокая точность** (поменяет `Real` на MPFR-50 и отключит каталог `tests`):
```bash
cmake -S . -B build -G Ninja -DUSE_MULTIPRECISION=ON
cmake --build build
```

**Тесты:**
```bash
ctest --test-dir build
```

**Примеры** собираются в `build/examples/`. Сборка использует флаги `-Wall -O3 -march=native`.

---

## Быстрый старт

Минимальный пример — маятник. Вектор состояния имеет размерность $2N$: первая половина — координаты $q$, вторая — сопряжённые импульсы $p$.

```cpp
#include <numbers>
#include <auto_diff.hpp>
#include <auto_diff_base.hpp>
#include <yoshida/leapfrog_extension.hpp>
#include <yoshida/yoshida.hpp>
#include <runge_kutta/odesolver.hpp>
using namespace AutoDiff;
using namespace boost::numeric::ublas;

int main() {
    auto q = Variable<0>{};
    auto p = Variable<1>{};
    constexpr double omega = 2 * std::numbers::pi;
    // H = p^2 / 2 - cos(omega * q)
    auto ham = p * p * Constant<0.5>{} - cos(Constant<omega>{} * q);

    // Leapfrog для расщеплённого гамильтониана + Yoshida 4-го порядка
    auto leapfrog   = Yoshida::make_leapfrog_for_sep<1>(ham);
    auto integrator = Yoshida::make_yoshida_from_leapfrog<4>(leapfrog);

    vector<Real> x(2);   // x[0] = q, x[1] = p
    x[0] = 0.25;
    x[1] = 0.0;

    RungeKutta::Integrate(x, 0.0, 2 * std::numbers::pi, 0.001, integrator);
}
```

Цепочка построения: **выражение гамильтониана** → `make_leapfrog_for_sep` (автоматическое расщепление) → `make_yoshida_from_leapfrog<order>` → `RungeKutta::Integrate`.

---

Ниже — справочник по отдельным модулям.

## precision

Параметризуемый тип вещественного числа и адаптеры для математических функций. Файл: `precision/precision.hpp`.

### Тип `Real`

По умолчанию `Real` — это обычный `double`:

```cpp
using Real = double;
```

При сборке с `-DUSE_MULTIPRECISION` тип переопределяется:

```cpp
using Real = boost::multiprecision::static_mpfr_float_50;
```

Это MPFR-тип с **фиксированной** точностью порядка 50 десятичных знаков — **не** произвольная точность; точность задаётся в момент компиляции.

## auto_diff

Модуль символической дифференциации, в котором математические выражения представляются **типами C++**, а дифференцирование выполняется на этапе компиляции: результат `make_diff` — это новое выражение (тип), а не численное значение. Файл-агрегатор — [`auto_diff/auto_diff.hpp`](auto_diff/auto_diff.hpp).

### Ядро: выражения как типы

Весь синтаксис строится на CRTP-базе `Expression<Derived>` (см. [`auto_diff/auto_diff_base.hpp`](auto_diff/auto_diff_base.hpp)). Дерево выражения — это дерево **типов**: `Add<A, B>`, `Mul<A, B>`, `Function<A, FuncTag>` и т.д. Базовый класс задаёт универсальный интерфейс:

```cpp
template <typename Derived> struct Expression {
    template <IsReal... Args> inline Real operator()(Args&&... x) const; // вычисление
    inline std::string print() const;        // строковое представление формулы
    inline size_t priority() const;          // приоритет для вставки скобок
};
```

Вызов `expr(x0, y0)` — это вычисление дерева в точке с аргументами типа `Real`.

### Базовые строительные блоки

* `Variable<Idx>` — переменная $x_{\text{Idx}}$; при вычислении берётся как `std::get<Idx>` из кортежа аргументов, печатается как `x_{Idx}`.
* `Constant<C>` — константа с `double`-значением `C`.
* `Integer<C>` — константа с целочисленным значением `C`.

```cpp
template <size_t Idx>  struct Variable : public Expression<Variable<Idx>>;
template <double Val>  struct Constant : public Expression<Constant<Val>>;
template <int Val>     struct Integer  : public Expression<Integer<Val>>;
```

### Операторы и функции

Операторы возвращают типы-узлы (см. [`auto_diff/auto_diff_expressions.hpp`](auto_diff/auto_diff_expressions.hpp)):

* `a + b` → `Add<A, B>`;
* `a - b` → `Add<A, Mul<Constant<-1.0>, B>>` (вычитание записано как сложение с умножением на $-1$); унарный минус `-a` → `Minus<Integer<0>, A>`;
* `a * b` → `Mul<A, B>`;
* `a / b` → `Mul<A, Div<Constant<1.0>, B>>` (деление — умножение на $1/b$).

Встроенные функции оборачивают sin, cos, exp и sqrt

```cpp
template <IsExpression Value> static Sin<Value>  sin(const Value& v);
template <IsExpression Value> static Cos<Value>  cos(const Value& v);
template <IsExpression Value> static Sqrt<Value> sqrt(const Value& v);
template <IsExpression Value> static Exp<Value>  exp(const Value& v);
```

Для произвольной пользовательской функции `Func` существует

```cpp
template <auto Func, IsExpression A>
constexpr static Function<A, FuntionTagFromFunc<Func>> apply_function(A a);
```

которая строит тип `Function<A, FuntionTagFromFunc<Func>>`.

### Концепты

Ключевые C++23-концепты (в [`auto_diff/auto_diff_base.hpp`](auto_diff/auto_diff_base.hpp)):

* `IsExpression<T>` — `T` наследуется от `Expression<T>`;
* `IsVariable<T>` — `T` есть `Variable<Index>`;
* `IsConstant<T>` — `T` есть `Constant<Value>`;
* `IsNumber<T>` — `IsInteger<T> || IsConstant<T>` (и `IsNotNumber<T>`).

Именно эти концепты определяют, в какие узлы разрешено собирать выражения.

### Дифференцирование: `make_diff`

```cpp
template <IsExpression A, IsVariable Var>
FixedPointSquise<typename Diff<A, Var>::type>::type make_diff(const A& a, const Var& b);
```

`Diff` — compile-time-трейт (см. [`auto_diff/auto_diff.hpp`](auto_diff/auto_diff.hpp)), реализующий стандартные правила дифференциирования через типы:

* $\dfrac{\partial x_{i}}{\partial x_{i}} = 1$, $\dfrac{\partial x_{i}}{\partial x_{j}} = 0$ при $i \neq j$; производные констант — $0$;
* сумма/разность — линейны;
* правило произведения: $\dfrac{\partial (AB)}{\partial x} = A\,\dfrac{\partial B}{\partial x} + \dfrac{\partial A}{\partial x}\,B$;
* правило частного: $\dfrac{\partial (A/B)}{\partial x} = \dfrac{\dfrac{\partial A}{\partial x}\,B - A\,\dfrac{\partial B}{\partial x}}{B^{2}}$;
* цепное правило для `Function<A, F>`: $\dfrac{\partial F(v)}{\partial x} = F'(v)\cdot \dfrac{\partial v}{\partial x}$.

Частные производные получаются простым выбором переменной (`make_diff(expr, Variable<0>{})`); производные высшего порядка — вложенным вызовом `make_diff(make_diff(f, x), y)`. Результат автоматически нормализован и упрощён (обёртка `FixedPointSquise`).

### Нормализация и упрощение

* `Normalize<T>` (в [`auto_diff/auto_diff_normalizer.hpp`](auto_diff/auto_diff_normalizer.hpp)) — разворачивает вложенные `Add`/`Mul` в плоские суммы/произведения и склеивает численные листья в одну константу (нейтральные элементы: $0$ для суммы, $1$ для произведения).
* `Squiser` — набор правил редукции (в [`auto_diff/auto_diff_squiser_rules.hpp`](auto_diff/auto_diff_squiser_rules.hpp)): $1\cdot A \to A$, $0 + A \to A$, $0\cdot A \to 0$, $A\cdot(1/A) \to 1$, $A - C \to -C + A$, $A + A \to 2A$, $0 - A \to -A$, $\sqrt{A}\cdot\sqrt{A} \to A$ и складывание численных констант. `Squise<T>` применяет правила рекурсивно ко всему дереву.
* `FixedPointSquise<T>` (в [`auto_diff/auto_diff_squise.hpp`](auto_diff/auto_diff_squise.hpp)) — итерация «Normalize + Squise» до неподвижной точки.
* `simplify` — runtime-обёртка, возвращающая упрощённый тип:

```cpp
template <IsExpression Expr>
typename FixedPointSquise<Expr>::type simplify(const Expr& expr);
```

### Подстановка

```cpp
template <IsExpression Expr, IsExpression Pattern, IsExpression Target>
constexpr typename Subst<Expr, Pattern, Target>::type
subst(const Expr& expr, const Pattern& pattern, const Target& target);
```

Заменяет все подвыражения, чей **тип** структурно равен `Pattern`, на `Target` (матчинг по типу, а не по значению). Функция `recursive_subst<N_iters>(expr, pattern, initial)` — итерация подстановки с упрощением на каждом шаге.

### `newton_expr`

Метод Ньютона, выполняемый в символьном пространстве:

```cpp
template <size_t N_iters, IsExpression Expr, IsExpression Init, size_t Idx>
constexpr decltype(auto) newton_expr(Expr expr, Variable<Idx> var, Init init);
```

Вычисляет $f^\prime = make\_diff(expr, var)$, строит символьную карту $x \mapsto x - \dfrac{f(x)}{f'(x)}$ и применяет её `N_iters` раз к начальному выражению `init`, упрощая результат после каждой итерации.

### Отладка: `print()`

`expr.print()` возвращает строковое представление формулы (например, `x_0 * x_1 + 2`); `priority()` управляет вставкой скобок при печати бинарных операций.

### Численные и аналитические производные

Встроенные аналитические производные заданы для `sin`, `cos`, `sqrt`, `exp` — это специализации шаблона `AutoDiff::AnalyticalDiff`. Для любой другой функции дифференцирование по умолчанию выполняется **численно** — центральной разностью при вычислении, с шагом $h = \sqrt{\varepsilon}$, где  $\varepsilon$ = `std::numeric_limits<Real>::epsilon()`:

$$
\frac{\partial F(v(x))}{\partial x} = \frac{F\big(v(x + h)\big) - F\big(v(x - h)\big)}{2h}
$$

Соответствующий тип `DefaultFuncDiff<A, Func, Idx>` хранит внутреннее выражение и вычисляет разность на лету.

### Аналитическая производная собственной функции

Механизм расширяется: аналитическую производную своей функции можно задать, определив **правило** — специализацию `AutoDiff::AnalyticalDiff`:

```cpp
Real my_func(Real x); // обычная функция 

namespace AutoDiff {
template <IsExpression A>
struct AnalyticalDiff<A, static_cast<Real (*)(Real)>(my_func)> {
    using type = /* выражение my_func'(A), напр. Mul<Constant<2.0>, A> */;
};
}
```

После этого `apply_function<my_func>(expr)` дифференцируется аналитически: выбирается перегрузка `Diff<Function<A, FuntionTagFromFunc<Func>>, Var>` с `AnalyticalDiff` (именно концепт `HasAnalyticalDiff` решает, какая ветка срабатывает), и производная строится по цепному правилу. Без специализации срабатывает численная ветка `DefaultFuncDiff`.

**Контракт по функции.** `Func` в `apply_function` передаётся как константное выражение и запекается в тип `FuntionTagFromFunc<Func>`. Поэтому в механизм проходят только **stateless** вызываемые объекты: обычные функции и лямбды **без захвата**. Лямбда с захватом или функтор со состоянием не подойдут — у них нет значения-константы. Это прямое следствие zero-cost природы `auto_diff`: вызываемый объект нигде не хранится в рантайме, он существует только как часть типа.

## Гамильтоновы системы (ham.hpp)

Модуль для построения правой части канонических уравнений из аналитического выражения гамильтониана. Файл: [`auto_diff/ham.hpp`](auto_diff/ham.hpp).

### `Hamiltonian::compile_system<N>`

```cpp
template <size_t N, typename H>
    requires IsExpression<H>
decltype(auto) compile_system(H ham);
```

Принимает аналитическое выражение гамильтониана `H` и возвращает замыкание

```cpp
(Real t, const vector<Real> &x, vector<Real> &result) -> void
```

которое вычисляет канонические уравнения:

$$
\dot{q} = \frac{\partial H}{\partial p}, \qquad \dot{p} = -\frac{\partial H}{\partial q},
$$

где производные берутся автоматически (через `make_diff`) по переменным `Variable<i+N>` (для $\dot{q}$) и `Variable<i>` (для $\dot{p}$).

### `call_with_vector_impl`

Вспомогательная функция, вычисляющая выражение на компонентах вектора состояния:

```cpp
template <typename Func, std::size_t... I>
decltype(auto) call_with_vector_impl(Real t,
                                      const boost::numeric::ublas::vector<Real> &vec,
                                      Func &&f, std::index_sequence<I...>);
```

Пакетом `vec[I]..., t` подставляет компоненты вектора в дифференцированное выражение.

### `IsRhs`

Возвращаемое замыкание удовлетворяет концепту `IsRhs`, поэтому его можно напрямую передать в конструкторы решателей:

```cpp
auto rhs    = Hamiltonian::compile_system<N>(ham);
auto solver = make_rk4(rhs);
auto brk    = make_brk(k, s, rhs);
```

### Вектор состояния

Размерность вектора состояния равна $2N$: первая половина содержит координаты $q_0, \dots, q_{N-1}$, вторая половина — сопряжённые импульсы $p_0, \dots, p_{N-1}$.

## newton

Численный блок для решения нелинейных систем и оптимизации. Все методы в файле [`newton/newton.hpp`](newton/newton.hpp).

### API

```cpp
int  solve(vector<Real> &x, std::function<void(const vector<Real> &, vector<Real> &)> f, Real damping = 1);
int  solve(Real &x, std::function<Real(Real)> f);
int  solve(Real &x, std::function<Real(Real)> f, std::function<Real(Real)> df);
int  solve_linear(const matrix<Real> &A, const vector<Real> &b, vector<Real> &x);
int  solve_linear(const matrix<Real> &A, const matrix<Real> &b, matrix<Real> &x);
int  inverse(const matrix<Real> &A, matrix<Real> &inverse);
Real find_minimum(Real a, Real b, std::function<Real(Real)> f);
```

* `solve` — метод Ньютона: векторная версия строит якобиан численно (разностно, с параллельным `#pragma omp parallel for`) и итеративно решает линейные системы; скалярная версия — с численной или аналитической производной.
* `solve_linear` — решение линейных систем с LU-разложением (`lu_factorize`/`lu_substitute`): правая часть — вектор или матрица.
* `inverse` — обращение матрицы через решение $A X = I$.
* `find_minimum` — одномерная минимизация методом золотой секции на отрезке $[a, b]$.

### Возврат

Все функции `solve*` и `inverse` возвращают `0` при успехе и ненулевое значение при неудаче (например, вырожденный якобиан или нулевая производная).

### Использование в проекте

В [`runge_kutta/include/runge_kutta/brk.hpp`](runge_kutta/include/runge_kutta/brk.hpp) модуль применяется для вычисления корней полиномов Лежандра — узлов Гаусса—Лежандра для схем интегрирования.

## runge_kutta

Модуль интегрирования ОДУ для гамильтоновых систем: набор Runge–Kutta методов (классический RK4 и семейство методов Батчера-Куцтнера ) с единым интерфейсом солвера и вспомогательными концептами. Все типы живут в неймспейсе `RungeKutta`.

### odesolver.hpp — общие контракты

Файл [`runge_kutta/include/runge_kutta/odesolver.hpp`](runge_kutta/include/runge_kutta/odesolver.hpp) задаёт три концепта:

* `IsRhs` — функция-«правая часть» системы $\dot{x} = f(t, x)$ со сигнатурой `(Real t, const vector<Real>& x, vector<Real>& result) -> void`, где `vector<Real>` — вектор из `boost::numeric::ublas`.
* `IODESolver` — любой солвер, реализующий `void Step(Real t, vector<Real>& x, Real h)`.
* `HasReset` — солвер, нуждающийся в сбросе внутреннего состояния через `void Reset()`.

На их основе определён свободный шаблон `Integrate` в двух перегрузках:

```cpp
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator);
void Integrate(vector<Real> &x, Real t0, Real t1, Real h, Solver &integrator,
               std::function<void(Real, const vector<Real> &)> callback);
```

Первая перегрузка просто выполняет шаги от $t_0$ до $t_1$ с шагом $h$. Вторая перед каждым шагом вызывает `callback(t, x)` и один раз — после последнего (с конечным состоянием). Обе вызывают `integrator.Reset()` перед циклом, если для `Solver` выполняется `HasReset` (через `if constexpr`), чтобы сбросить состояние, накопленное при предыдущей интеграции.

### rk4.hpp — классический RK4

Файл [`runge_kutta/include/runge_kutta/rk4.hpp`](runge_kutta/include/runge_kutta/rk4.hpp) содержит `RK4<Rhs>` и фабрику `make_rk4`:

```cpp
template <typename Rhs>
    requires IsRhs<Rhs>
RK4(int n, Rhs &&rhs);

template <typename Rhs>
    requires IsRhs<Rhs>
auto make_rk4(int n, Rhs &&rhs);
```

`n` — размерность системы (рантаймовый параметр; размер векторов `k1..k4`, выделяемых при конструировании). Метод `Step` — стандартная четырёхстадийная схема:

$$
x \leftarrow x + \frac{h}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right),
$$

где $k_i$ — оценки производной в узлах $t$, $t + h/2$ и $t + h$.

### brk.hpp — Методы Батчера-Куцтнера

Файл [`runge_kutta/include/runge_kutta/brk.hpp`](runge_kutta/include/runge_kutta/brk.hpp) содержит `BRK<k, s, Rhs>` и фабрику `make_brk`:

```cpp
template <size_t k, size_t s, typename Rhs>
    requires IsRhs<Rhs>
BRK(int n, Rhs &&rhs);

template <size_t k, size_t s, typename Rhs>
    requires IsRhs<Rhs>
auto make_brk(int n, Rhs &&rhs);   // вызов: make_brk<2, 3>(n, rhs)
```

Здесь $k$ и $s$ — шаблонные параметры: $k$ — число узлов сетки Гаусса-Лежандра, $s$ — число внутренних подшагов метода простых итераций на каждом шаге интегрирования. Размерность системы $n$. Узлы $c_i$ и матрицы переноса $A$, $b$, $L$ считаются в рантайме: корни многочлена Лежандра $L_k(x)$ находятся методом Ньютона (`Newton::solve`), веса и матрицы строятся из этих корней.

Схема шага: по $k$ узлам оценивается правая часть $\dot{x}$ в точках $t + c_i h$, затем $s$ раз применяются переносы по $A$, и состояние обновляется как $x \leftarrow x + h\, K b$. Первый шаг особый: его подготовка делается внутренним вызовом RK4 — метод `Prepare(t, x, h)` шагает RK4 от $t$ до каждого узла $t + c_i h$. Флаг `_first_step` (сбрасывается `Reset()`) включает этот режим ровно на первом шаге; `Integrate` вызывает `Reset()` автоматически, если солвер поддерживает `HasReset`.

### Когда что использовать

* **RK4** — универсальный метод 4-го порядка; достаточен для большинства задач с умеренными требованиями к точности, не требует подготовки узлов.
* **BRK** — симплектический метод высокого порядка (порядка $2k-1$ при $k$ узлах); применяется для решения гамильтоновых систем общего вида. Узлы $c_i$ и матрицы $A$, $b$, $L$ считаются в конструкторе методом Ньютона.

## parallel

Обёртки над OpenMP для параллельных циклов и редуцирования. Файл: [`parallel/parallel.hpp`](parallel/parallel.hpp).

### Концепты

```cpp
template <typename U, typename T>
concept IndexableContainer = requires(U u, size_t i) {
    { u[i] } -> std::same_as<T &>;
    { u.size() } -> std::same_as<size_t>;
};

template <typename Op, typename T>
concept BinaryOperation = requires(Op &&op, const T &a, const T &b) {
    { op(a, b) } -> std::same_as<T>;
};
```

`IndexableContainer` требует доступа по индексам `u[i] -> T&` и `u.size() -> size_t`; `BinaryOperation` — операцию `op(a, b) -> T`.

### `for_i`

Параллельный цикл по диапазону индексов:

```cpp
template <typename Fn>
    requires std::invocable<Fn &, size_t>
void for_i(Fn &&f, int start, int end);
```

### `for_each`

Параллельный проход по элементам коллекции:

```cpp
template <typename T, typename Fn, typename U>
    requires std::invocable<Fn &, const T &> && IndexableContainer<U, T>
void for_each(Fn &&f, U &collection);
```

### `reduce`

Параллельное сворачивание (fold) по коллекции ассоциативной бинарной операцией с заданным нейтральным элементом:

```cpp
template <typename T, typename Op, typename U>
    requires BinaryOperation<Op, T> && IndexableContainer<U, T>
T reduce(Op &&binary_op, const U &collection, T identity);
```

### Требования

Модуль требует OpenMP. Число потоков везде ограничено `omp_get_max_threads()` и размером работающей области — параллелизм включается только когда объёма работы достаточно.

## Yoshida

Семейство методов Йошиды и метод Leapfrog рассчитаны на расщепленные гамильтонианы вида:
$$H(p, q, t) = H_1(p, q, t) + H_2(p, q, t)$$

Можно выделить 4 специальных типа расщеплений:
* **I.** $H(p, q, t) = H_1(p) + H_2(q, t);$
* **II.** $H(p, q, t) = H_1(p, q) + H_2(q, t),$ где $\cfrac{\partial^2 H_1}{\partial p \partial q} \not\equiv 0$ и $H_1$ интегрируем;
* **III.** $H(p, q, t) = H_1(p, q) + H_2(p, t),$ где $\cfrac{\partial^2 H_1}{\partial p \partial q} \not\equiv 0$ и $H_1$ интегрируем;
* **IV.** $H(p, q) = H_1(p, q) + H_2(p, q),$ где $\cfrac{\partial^2 H_1}{\partial p \partial q} \not\equiv 0$, $\cfrac{\partial^2 H_2}{\partial p \partial q} \not\equiv 0$ и $H_1, H_2$ интегрируемы.

Методы Йошиды порождаются из соответствующего Leapfrog с помощью функции:
```cpp
template <size_t order, typename h1_solver, typename h2_solver>
auto make_yoshida_from_leapfrog(LeapFrog<h1_solver, h2_solver> leapfrog); 
```

Для создания метода Leapfrog представлены четыре функции под различные типы гамильтонианов.

### Общие требования к солверам (C++20 Concept)

Некоторые фабричные функции принимают пользовательские функторы (`h1`, `h2`), которые должны предоставлять точные решения подсистем и удовлетворять концепту `IsSolverFunctor`:

```cpp
template <typename F>
concept IsSolverFunctor = requires(F f, Real t, boost::numeric::ublas::vector<Real> &x, Real h) {
    { f(t, x, h) } -> std::same_as<void>;
};
```

#### Параметры сигнатуры функтора:
* **`t`** (`Real`) — текущая временная метка (time step). Передается по значению.
* **`x`** (`boost::numeric::ublas::vector<Real>&`) — неконстантная ссылка на вектор состояния. Имеет четную размерность: первую половину составляют обобщенные координаты ($q$), а вторую — сопряженные к ним импульсы ($p$). Функтор модифицирует этот вектор напрямую (in-place).
* **`h`** (`Real`) — величина шага интегрирования (step size). Функтор должен вернуть решение в точке $t + h$, используя начальные данные из точки $t$.
* **Возвращаемое значение:** Должно быть **строго** `void`.

### Фабричные функции Leapfrog

#### 1. Для гамильтонианов I типа (автоматическое расщепление)
Определена в заголовочном файле [`yoshida/include/yoshida/leapfrog_extension.hpp`](yoshida/include/yoshida/leapfrog_extension.hpp):
```cpp
template <size_t NFreedom, typename H>
    requires IsExpression<H>
decltype(auto) make_leapfrog_for_sep(const H &ham);
```

#### 2. Для гамильтонианов II типа
Определена в заголовочном файле [`yoshida/include/yoshida/leapfrog_extension.hpp`](yoshida/include/yoshida/leapfrog_extension.hpp):
```cpp
template <size_t NFreedom, typename H, typename h2_solver>
    requires IsExpression<H> && IsSolverFunctor<h2_solver>
decltype(auto) make_leapfrog_for_spv(const H &ham, h2_solver &&h2);
```
* **`h2`** — функтор, предоставляющий решение для «нестандартной» (не расщепляемой автоматически) части гамильтониана ($H_2$); вторая часть Leapfrog строится автоматически из `ham`. `ham` это $H_1$.

#### 3. Для гамильтонианов III типа
Определена в заголовочном файле [`yoshida/include/yoshida/leapfrog_extension.hpp`](yoshida/include/yoshida/leapfrog_extension.hpp):
```cpp
template <size_t NFreedom, typename H, typename h2_solver>
    requires IsExpression<H> && IsSolverFunctor<h2_solver>
decltype(auto) make_leapfrog_for_spq(const H &ham, h2_solver &&h2);
```
* **`h2`** — функтор, предоставляющий решение для «нестандартной» части гамильтониана ($H_2$); другая часть Leapfrog строится автоматически из `ham`. `ham` это $H_1$.

#### 4. Универсальный интерфейс (для гамильтонианов IV типа)
Определена в заголовочном файле [`yoshida/include/yoshida/leapfrog.hpp`](yoshida/include/yoshida/leapfrog.hpp):
```cpp
template <typename h1_solver, typename h2_solver>
    requires IsSolverFunctor<h1_solver> && IsSolverFunctor<h2_solver>
auto make_leapfrog(h1_solver &&h1, h2_solver &&h2);
```
* **`h1`** — функтор, предоставляющий аналитическое или численное решение для гамильтониана $H_1$.
* **`h2`** — функтор, предоставляющий решение для гамильтониана $H_2$.

### Порядок метода

* Метод Йошиды порождаётся **рекурсивно** из Leapfrog: базовый случай — порядок $2$ (чистый Leapfrog), каждый следующий порядок складывается из трёх шагов «вложенного» метода порядка на два меньше. Поэтому поддерживаются **только чётные порядки**: $2, 4, 6, 8, \dots$ (например, `make_yoshida_from_leapfrog<4>` — схема 4-го порядка).
* Коэффициенты схемы `YoshidaCoeffs<order>` вычисляются один раз итеративно (метод Ньютона, 10 итераций) в рантайме.

## Примеры

### pendulum

В примере [`examples/pendulum/pendulum.cpp`](examples/pendulum/pendulum.cpp) рассматривается маятник с частотой $\omega = 2\pi$ и гамильтонианом

$$
H = \frac{p^2}{2} - \cos(2\pi q).
$$

Гамильтониан расщеплён на кинетическую $H_1 = p^2/2$ и потенциальную $H_2 = -\cos(2\pi q)$ части, поэтому интегратор строится из `make_leapfrog_for_sep` (Leapfrog для расщеплённого гамильтониана I типа) и `make_yoshida_from_leapfrog<4>` (схема Йошиды 4-го порядка):

```cpp
auto ham = p * p * Constant<0.5>{} - cos(Constant<omega>{} * q);
auto leapfrog = Yoshida::make_leapfrog_for_sep<1>(ham);
auto integrator = Yoshida::make_yoshida_from_leapfrog<4>(leapfrog);
```

Партия из $N = 20$ траекторий с различными начальными условиями интегрируется за один период $2\pi$ параллельно через `Parallel::for_i`; результаты пишутся в `traj.dat`. Файл `plot_scr` — скрипт gnuplot, который строит фазовые траектории на плоскости $(q, p)$ и проверяет сохранение энергии: каждая кривая подписывается значением $H$ вдоль траектории.

Запуск: `./build/examples/pendulum/pendulum`, затем `gnuplot plot_scr` (скрипт копируется рядом с бинарником при сборке).

### sitnikov

Пример [`examples/sitnikov/sitnikov.cpp`](examples/sitnikov/sitnikov.cpp) — построение карты отображения Пуанкаре для задачи Ситникова. Гамильтониан явно зависит от времени:

$$
H(p, q, t) = \frac{p^2}{2} - \frac{1}{\sqrt{r^2(t) + q^2}},
$$

где $q$ — координата третьего тела, $p$ — сопряжённый импульс, $r^2(t)$ — квадрат расстояния между двумя фиксированными телами. Время $t$ передается как дополнительная координата — `Variable<2>` — а пользовательская функция вставляется в выражение через `apply_function`:

```cpp
auto t = Variable<2>{};
auto r = apply_function<distance_sqr_ch>(t);
auto z = sqrt(r + q * q);
auto ham = p * p * Constant<0.5>{} - Constant<1.0>{} / z;
```

`distance_sqr_ch` — кэшированная версия `distance_sqr`: расстояние предварительно вычисляется (решение уравнения Кеплера) в массив из 10000 точек за период и далее оценивается линейной интерполяцией, что делает вычисление $H$ в каждом шаге интегрирования быстрым.

Пример демонстрирует поддержку явно зависящих от времени гамильтонианов $H(p, q, t)$, встраивание пользовательских функций в выражения и крупный параллельный пакет: начальными условиями заполняется сетка $128 \times 128$, каждая траектория интегрируется за период $2\pi$ (шаг $0.1$) параллельно через `Parallel::for_i`, а результаты пишутся в `out.dat`.

### three_body

В примере [`examples/three_body/three_body.cpp`](examples/three_body/three_body.cpp) задача $N$ тел собирается с помощью заголовка [`examples/three_body/nbody.hpp`](examples/three_body/nbody.hpp):

* `Body<mass>` — тело, масса которого задана шаблонным параметром;
* `make_planar_nbody` / `make_general_nbody` — сборка задачи в 2D/3D по списку тел;
* `NBodyProblem::compiled_system()` — скомпилированная система дифференциальных уравнений (через `Hamiltonian::compile_system<NFreedom>`).

Гамильтониан — сумма кинетической энергии и Ньютоновского потенциала:

$$
H = \sum_i \frac{|\mathbf{p}_i|^2}{2 m_i} - \sum_{i < j} \frac{1}{|\mathbf{r}_i - \mathbf{r}_j|}.
$$

Интегрирование выполняется методом Рунге—Кутты из `runge_kutta/brk.hpp`:

```cpp
auto body1 = Body<1.0>{};
auto body2 = Body<1.0>{};
auto body3 = Body<1.0>{};
auto problem = make_planar_nbody(body1, body2, body3);
auto system = problem.compiled_system();
auto integrator = RungeKutta::make_brk<2, 3>(problem.NCoords, system);
```

Запуск: `./build/examples/three_body/three_body`.

## Тесты

Тесты написаны на **doctest** и регистрируются в **CTest** через `doctest_discover_tests` (см. [`tests/CMakeLists.txt`](tests/CMakeLists.txt)). Все тесты находятся в каталоге `tests/auto_diff/`:

* **`basic_rules`** — базовые правила дифференцирования: производная константы, переменной, правила суммы, произведения и частного.
* **`chain_rule`** — правило цепи для вложенных составных функций.
* **`partial`** — частные производные по нескольким переменным и равенство смешанных вторых частных производных (теорема Шварца).
* **`random_expr`** — для случайно сгенерированных выражений (генерация в момент компиляции по сидам) символьная производная сверяется с центральной конечной разностью.

Запуск:

```
ctest --test-dir build
```

Примечание: тесты собираются только при отсутствии определения `USE_MULTIPRECISION` (каталог `tests` добавляется в [`CMakeLists.txt`](CMakeLists.txt) под условием `if(NOT DEFINED USE_MULTIPRECISION)`).

## Лицензия и стиль

* **Лицензия** — Mozilla Public License 2.0 (см. [`LICENSE`](LICENSE)).
* **Стиль кода** — см. [`.clang-format`](.clang-format): база LLVM, отступ 4 пробела, без табуляции, `NamespaceIndentation: Inner`.
* Сборка компилируется с `-Wall`; в пользовательском коде предпочтительно использовать `Real` (см. [precision](#precision)).
