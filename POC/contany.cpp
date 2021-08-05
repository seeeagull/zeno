#include <iostream>
#include <variant>
#include <cstdint>
#include <any>

using scalar_type_variant = std::variant
        < bool
        , uint8_t
        , uint16_t
        , uint32_t
        , uint64_t
        , int8_t
        , int16_t
        , int32_t
        , int64_t
        , float
        , double
        >;

template <class T, class = void>
struct any_traits {
    using underlying_type = T;
};

template <class T>
struct any_traits<T, std::void_t<decltype(
        std::declval<scalar_type_variant &>() = std::declval<T>()
        )>> {
    using underlying_type = scalar_type_variant;
};

struct any : std::any {
    any() = default;
    any(any &&a) = default;
    any(any const &a) = default;

    template <class T>
    any(T const &t)
    : std::any(static_cast<typename any_traits<T>::underlying_type>(t))
    {
        printf("!!!%s\n", typeid(t).name());
    }

    template <class T>
    any &operator=(T const &t) {
        std::any::operator=(
                static_cast<typename any_traits<T>::underlying_type>(t));
        return *this;
    }
};

template <class T>
struct is_variant : std::false_type {
};

template <class ...Ts>
struct is_variant<std::variant<Ts...>> : std::true_type {
};

template <class T>
auto our_any_cast(any const &a) {
    auto &&v = std::any_cast<typename any_traits<T>::underlying_type>(a);
    if constexpr (is_variant<std::decay_t<decltype(v)>>::value) {
        return std::get<T>(v);
    } else {
        return v;
    }
}

int main() {
    any a = 3.14f;
    std::cout << our_any_cast<int>(a) << std::endl;
    return 0;
}
