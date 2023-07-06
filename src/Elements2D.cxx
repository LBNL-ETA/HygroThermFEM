#include <numeric>

#ifdef STL_MULTITHREADING
#    include <execution>
#endif

#include "Elements2D.hxx"
#include "NodePool.hxx"

namespace HygroThermFEM
{
    SquareMatrix ElementsLinear2D::conductanceMatrix()
    {
        SquareMatrix result{NodePool::Instance().maxIndex()};

#ifdef STL_MULTITHREADING

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
                                  result(indexes[i] - 1, indexes[j] - 1) +=
                                    conductance(i, j) + condDer(i, j);
                              }
                          }
                          mtx.unlock();
                      });

#else
        for(const auto & element : m_Elements)
        {
            auto indexes = element->nodeIndexes();
            auto conductance = element->DDuMatrices();
            // auto testConductance = conductance.toVector();
            auto condDer = element->DpDuMatrices();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                {
                    result(indexes[i] - 1, indexes[j] - 1) += conductance(i, j) + condDer(i, j);
                }
            }
        }
#endif

        return result;
    }

    std::vector<double> ElementsLinear2D::getLumpedMass(const double DTime)
    {
        const auto numOfNodes{NodePool::Instance().maxIndex()};
        std::vector<std::vector<double>> Capacitance(numOfNodes,
                                                     std::vector<double>(numOfNodes, 0));

#ifdef STL_MULTITHREADING
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
#else
        for(const auto & element : m_Elements)
        {
            auto indexes = element->nodeIndexes();
            auto capacitance = element->capacitanceMatrices();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                {
                    Capacitance[indexes[i] - 1][indexes[j] - 1] += capacitance(i, j);
                }
            }
        }
#endif

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

#ifdef STL_MULTITHREADING

        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto capacitance = aElement->capacitanceMatrices();
                          mtx.lock();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  Capacitance(indexes[i] - 1, indexes[j] - 1) +=
                                    capacitance(i, j) / DTime;
                              }
                          }
                          mtx.unlock();
                      });

#else
        for(const auto & element : m_Elements)
        {
            auto indexes = element->nodeIndexes();
            auto capacitance = element->capacitanceMatrices();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                {
                    Capacitance(indexes[i] - 1, indexes[j] - 1) += capacitance(i, j) / DTime;
                }
            }
        }
#endif

        return SquareMatrix{Capacitance};
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

#ifdef STL_MULTITHREADING

        std::mutex mtx;

        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          const auto indexes = aElement->nodeIndexes();
                          const auto vecR = aElement->rightSideVector();
                          mtx.lock();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              result[indexes[i] - 1] += vecR[i];
                          }
                          mtx.unlock();
                      });

#else
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto vecR = element->rightSideVector();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                result[indexes[i] - 1] += vecR[i];
            }
        }
#endif

        return result;
    }

    std::vector<NodeFlux> ElementsLinear2D::flux() const
    {
        std::vector<NodeFlux> result(NodePool::Instance().maxIndex(), {0, 0});

        std::vector<std::vector<NodeFlux>> fluxes(NodePool::Instance().maxIndex(),
                                                  std::vector<NodeFlux>());

#ifdef STL_MULTITHREADING
        std::mutex mtx;

        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          const auto indexes = aElement->nodeIndexes();
                          const auto flux = aElement->flux();
                          mtx.lock();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              fluxes[indexes[i] - 1].push_back(flux[i]);
                          }
                          mtx.unlock();
                      });

#else
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto flux = element->flux();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                fluxes[indexes[i] - 1].push_back(flux[i]);
            }
        }
#endif

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

    void ElementsLinear2D::clearElements()
    {
        m_Elements.clear();
    }


}   // namespace HygroThermFEM
