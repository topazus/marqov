#include <iostream>
#include <type_traits>
#include <tuple>
#include <functional>

using namespace std;

struct Duck {};

template <typename... Ts>
struct Test : std::false_type
{};

template <typename Lattice, class... Ts>
struct Test<std::tuple<Lattice&, Ts...> , Lattice> : std::true_type
{
};

int main()
{
std::string t("test");
Duck duck;
Duck& duckref = duck;
auto inp = tuple_cat(std::forward_as_tuple(duckref), make_tuple(5, 5, 5.0));
//int j = inp
// int& i = typename std::tuple_element<0, decltype(inp)>::type();
std::cout<<std::boolalpha<< typename Head_is_ref_of_duck<decltype(inp), Duck>::Result ()<<std::endl;
std::cout<<Test<decltype(inp), string>()<<std::endl;
}