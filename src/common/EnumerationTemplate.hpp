#pragma once

namespace HygroThermFEM
{
    template<typename T>
    class Enum
    {
    public:
        class Iterator
        {
        public:
            explicit Iterator(const int value) : m_value(value)
            {}

            T operator*()const
            {
                return static_cast<T>(m_value);
            }

            void operator++()
            {
                ++m_value;
            }

            bool operator!=(Iterator rhs)
            {
                return m_value != rhs.m_value;
            }

        private:
            int m_value;
        };
    };
}