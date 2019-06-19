#include <numeric>
#include <execution>

#include "Elements2D.hxx"
#include "NodePool.hxx"

namespace HygroThermFEM
{
    SquareMatrix ElementsLinear2D::conductanceMatrix()
    {
        const auto numOfNodes{NodePool::Instance().maxIndex()};
        std::vector<std::vector<double>> result{numOfNodes, std::vector<double>(numOfNodes, 0)};
        // SquareMatrix result{NodePool::Instance().maxIndex()};
        // now integrate element matrices into global matrix

        std::mutex mtx;

        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto conductance = aElement->DDuMatrices();
                          // auto testConductance = conductance.toVector();
                          auto condDer = aElement->DpDuMatrices();
                          // auto testCondDer = condDer.toVector();
                          mtx.lock();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  result[indexes[i] - 1][indexes[j] - 1] +=
                                    conductance(i, j) + condDer(i, j);
                              }
                          }
                          mtx.unlock();
                      });
        return SquareMatrix{result};
    }

    std::vector<double> ElementsLinear2D::getLumpedMass(const double DTime)
    {
        const auto numOfNodes{NodePool::Instance().maxIndex()};
        std::vector<std::vector<double>> Capacitance(numOfNodes,
                                                     std::vector<double>(numOfNodes, 0));

        std::mutex mtx;

        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto capacitance = aElement->capacitanceMatrices();
                          // auto capTest = capacitance.toVector();
                          mtx.lock();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  Capacitance[indexes[i] - 1][indexes[j] - 1] += capacitance(i, j);
                              }
                          }
                          mtx.unlock();
                      });

        const auto size = Capacitance.size();

        std::vector<double> M(size, 0);

        // Creates lump matrix
        for(size_t i = 0; i < size; ++i)
        {
            for(size_t j = 0; j < size; ++j)
            {
                M[i] += Capacitance[i][j];
            }
            M[i] /= DTime;
        }

        return M;
    }

    SquareMatrix ElementsLinear2D::getMassMatrix(const double DTime)
    {
        SquareMatrix Capacitance{NodePool::Instance().maxIndex()};

        // now integrate element matrices into global matrix
        for(auto & aElement : m_Elements)
        {
            auto indexes = aElement->nodeIndexes();
            auto capacitance = aElement->capacitanceMatrices();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                {
                    Capacitance(indexes[i] - 1, indexes[j] - 1) += capacitance(i, j) / DTime;
                }
            }
        }

        return Capacitance;
    }

    bool ElementsLinear2D::isLinear() const
    {
        bool isLinear = true;
        for(auto & elem : m_Elements)
        {
            isLinear = isLinear && elem->isLinear();
            if(!isLinear)   // no need to waste time in loop
            {
                break;
            }
        }
        return isLinear;
    }

    void ElementsLinear2D::updateNodeValues(const std::vector<double> & values,
                                            const BaseVariable property,
                                            bool updatePreviousValue)
    {
        for(auto & aElement : m_Elements)
        {
            for(auto i = 0u; i < numOfQuadrilateralNodes; ++i)
            {
                auto & node = aElement->getNode(i);
                const auto index = node.getNodeNumber();
                node.setStateProperty(property, values[index - 1], updatePreviousValue);
            }
        }
    }

    IElementLinear2D * ElementsLinear2D::findElement(const size_t index1, const size_t index2)
    {
        IElementLinear2D * el = nullptr;
        for(auto & element : m_Elements)
        {
            if(element->haveBothNodes(index1, index2))
            {
                el = element.get();
            }
        }
        return el;
    }

    void ElementsLinear2D::assignElement(std::unique_ptr<IElementLinear2D> && el)
    {
        m_Elements.push_back(std::move(el));
    }

    std::vector<double> ElementsLinear2D::RVector() const
    {
        std::vector<double> result(NodePool::Instance().maxIndex(), 0);
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto vecR = element->rightSideVector();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                result[indexes[i] - 1] += vecR[i];
            }
        }
        return result;
    }

    std::vector<NodeFlux> ElementsLinear2D::flux() const
    {
        std::vector<NodeFlux> result(NodePool::Instance().maxIndex(), {0, 0});

        std::vector<std::vector<NodeFlux>> fluxes(NodePool::Instance().maxIndex(),
                                                  std::vector<NodeFlux>());

        // First pickup all fluxes from elements
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto flux = element->flux();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                fluxes[indexes[i] - 1].push_back(flux[i]);
            }
        }

        // Now need to average them
        for(size_t j = 0; j < fluxes.size(); ++j)
        {
            const double x = std::accumulate(fluxes[j].begin(),
                                             fluxes[j].end(),
                                             0.0,
                                             [&](double lhs, NodeFlux & a) { return lhs + a.x; })
                             / fluxes[j].size();
            const double y = std::accumulate(fluxes[j].begin(),
                                             fluxes[j].end(),
                                             0.0,
                                             [&](double lhs, NodeFlux & a) { return lhs + a.y; })
                             / fluxes[j].size();
            result[j] = {x, y};
        }

        return result;
    }

    const std::vector<std::unique_ptr<IElementLinear2D>> & ElementsLinear2D::elements() const
    {
        return m_Elements;
    }


}   // namespace HygroThermFEM
