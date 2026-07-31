#pragma once

#include <array>
#include <cstddef>

namespace HygroThermFEM
{
    //! \brief Small inline container for the per-node values of one element or boundary segment.
    //!
    //! An element carries four nodes and a boundary segment two, so the capacity is fixed and the
    //! storage can live inline. This matters because a container of this shape is produced once
    //! per coefficient term per element on every matrix assembly, and assemblies run several times
    //! per Newton-Raphson iteration: returning a std::vector meant a heap allocation and release
    //! for what is at most four doubles.
    template<typename ValueType>
    class NodalArray
    {
    public:
        static constexpr std::size_t maxNodes{4u};

        NodalArray() = default;

        //! Sized construction. Entries are value-initialised.
        explicit NodalArray(const std::size_t count) : m_Count(count)
        {}

        [[nodiscard]] std::size_t size() const
        {
            return m_Count;
        }

        [[nodiscard]] bool empty() const
        {
            return m_Count == 0u;
        }

        ValueType & operator[](const std::size_t index)
        {
            return m_Values[index];
        }

        const ValueType & operator[](const std::size_t index) const
        {
            return m_Values[index];
        }

        [[nodiscard]] ValueType * data()
        {
            return m_Values.data();
        }

        [[nodiscard]] const ValueType * data() const
        {
            return m_Values.data();
        }

        [[nodiscard]] ValueType * begin()
        {
            return m_Values.data();
        }

        [[nodiscard]] ValueType * end()
        {
            return m_Values.data() + m_Count;
        }

        [[nodiscard]] const ValueType * begin() const
        {
            return m_Values.data();
        }

        [[nodiscard]] const ValueType * end() const
        {
            return m_Values.data() + m_Count;
        }

        void push_back(const ValueType & value)
        {
            m_Values[m_Count] = value;
            ++m_Count;
        }

    private:
        std::array<ValueType, maxNodes> m_Values{};
        std::size_t m_Count{0u};
    };

    //! Coefficient values at the nodes of one element or boundary segment.
    using NodalValues = NodalArray<double>;
}   // namespace HygroThermFEM
