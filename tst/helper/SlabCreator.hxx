#pragma once

#include <vector>
#include <string>

#include "HygroThermFEM2D.hxx"

namespace TestHelper
{
    //! Starting corner for element node ordering
    //! Node positions for element at index i:
    //!   BottomLeft=2i-1, TopLeft=2i, BottomRight=2i+1, TopRight=2i+2
    enum class StartCorner
    {
        BottomLeft,
        BottomRight,
        TopLeft,
        TopRight
    };

    //! Direction for traversing element nodes
    enum class Direction
    {
        Clockwise,
        CounterClockwise
    };

    //! \brief Builder for creating simple slab geometries in tests
    //!
    //! Usage:
    //!   SlabBuilder(multiDomain)
    //!       .gridXCoordinates({0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0})
    //!       .height(0.05)
    //!       .material(material.name())
    //!       .state(state)
    //!       .build();
    class SlabBuilder
    {
    public:
        explicit SlabBuilder(HygroThermFEM::MultiDomain & domain) : m_Domain(domain)
        {}

        SlabBuilder & gridXCoordinates(std::vector<double> coords)
        {
            m_GridXCoordinates = std::move(coords);
            return *this;
        }

        SlabBuilder & height(double val)
        {
            m_Height = val;
            return *this;
        }

        SlabBuilder & material(std::string name)
        {
            m_Material = std::move(name);
            return *this;
        }

        SlabBuilder & state(HygroThermFEM::State val)
        {
            m_State = val;
            return *this;
        }

        SlabBuilder & startCorner(StartCorner corner)
        {
            m_StartCorner = corner;
            return *this;
        }

        SlabBuilder & direction(Direction dir)
        {
            m_Direction = dir;
            return *this;
        }

        void build()
        {
            // Create nodes
            for(const auto & xVal : m_GridXCoordinates)
            {
                m_Domain.nodes().createNode({.x = xVal, .y = 0.0, .state = m_State});
                m_Domain.nodes().createNode({.x = xVal, .y = m_Height, .state = m_State});
            }

            // Create elements
            const auto numElements = (m_Domain.nodes().maxIndex() - 2) / 2;
            for(size_t idx = 1; idx <= numElements; ++idx)
            {
                const size_t botLeft = 2 * idx - 1;
                const size_t topLeft = 2 * idx;
                const size_t botRight = 2 * idx + 1;
                const size_t topRight = 2 * idx + 2;

                // Order nodes based on start corner and direction
                // CCW from corner: corner → next CCW → opposite → prev CCW
                // CW from corner: corner → next CW → opposite → prev CW
                size_t node1{0}, node2{0}, node3{0}, node4{0};
                const bool isCW = (m_Direction == Direction::Clockwise);

                switch(m_StartCorner)
                {
                    case StartCorner::BottomLeft:
                        node1 = botLeft;
                        node2 = isCW ? topLeft : botRight;
                        node3 = topRight;
                        node4 = isCW ? botRight : topLeft;
                        break;
                    case StartCorner::BottomRight:
                        node1 = botRight;
                        node2 = isCW ? botLeft : topRight;
                        node3 = topLeft;
                        node4 = isCW ? topRight : botLeft;
                        break;
                    case StartCorner::TopLeft:
                        node1 = topLeft;
                        node2 = isCW ? topRight : botLeft;
                        node3 = botRight;
                        node4 = isCW ? botLeft : topRight;
                        break;
                    case StartCorner::TopRight:
                        node1 = topRight;
                        node2 = isCW ? botRight : topLeft;
                        node3 = botLeft;
                        node4 = isCW ? topLeft : botRight;
                        break;
                }

                m_Domain.createElement({.node1 = node1,
                                        .node2 = node2,
                                        .node3 = node3,
                                        .node4 = node4,
                                        .material = m_Material});
            }
        }

    private:
        HygroThermFEM::MultiDomain & m_Domain;
        std::vector<double> m_GridXCoordinates;
        double m_Height = 0.05;
        std::string m_Material;
        HygroThermFEM::State m_State = {};
        StartCorner m_StartCorner = StartCorner::BottomLeft;
        Direction m_Direction = Direction::CounterClockwise;
    };

}   // namespace TestHelper
