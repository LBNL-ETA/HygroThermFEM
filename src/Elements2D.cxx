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
        SquareMatrix result{HygroThermFEM::maxNodeIndex()};

#ifdef STL_MULTITHREADING

        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto conductance = aElement->DDuMatrices();
                          auto condDer = aElement->DpDuMatrices();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  std::lock_guard guard(mtx);
                                  result(indexes[i] - 1, indexes[j] - 1) +=
                                    conductance(i, j) + condDer(i, j);
                              }
                          }
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
        const auto numOfNodes{HygroThermFEM::maxNodeIndex()};
        std::vector<std::vector<double>> Capacitance(numOfNodes,
                                                     std::vector<double>(numOfNodes, 0));

#ifdef STL_MULTITHREADING
        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto capacitance = aElement->capacitanceMatrices();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  std::lock_guard guard(mtx);
                                  Capacitance[indexes[i] - 1][indexes[j] - 1] += capacitance(i, j);
                              }
                          }
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

    SquareMatrix ElementsLinear2D::getCMatrix()
    {
        SquareMatrix Capacitance{HygroThermFEM::maxNodeIndex()};

#ifdef STL_MULTITHREADING

        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto capacitance = aElement->capacitanceMatrices();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  std::lock_guard guard(mtx);
                                  Capacitance(indexes[i] - 1, indexes[j] - 1) +=
                                    capacitance(i, j);
                              }
                          }
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
                    Capacitance(indexes[i] - 1, indexes[j] - 1) += capacitance(i, j);
                }
            }
        }
#endif

        // Creates lump matrix
        SquareMatrix Mlump(HygroThermFEM::maxNodeIndex());
        Mlump.setZeros();
        for(size_t i = 0; i < HygroThermFEM::maxNodeIndex(); ++i)
        {
            for(size_t j = 0; j < HygroThermFEM::maxNodeIndex(); ++j)
            {
                Mlump(i,i) += Capacitance(i,j);
            }
        }

        return Mlump;
    }

    SquareMatrix ElementsLinear2D::getMassMatrix(const double DTime)
    {
        SquareMatrix Capacitance{HygroThermFEM::maxNodeIndex()};

#ifdef STL_MULTITHREADING

        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          auto indexes = aElement->nodeIndexes();
                          auto capacitance = aElement->capacitanceMatrices();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  std::lock_guard guard(mtx);
                                  Capacitance(indexes[i] - 1, indexes[j] - 1) +=
                                    capacitance(i, j) / DTime;
                              }
                          }
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
        for(auto & elem : m_Elements)
        {
            if(!elem->isLinear())
            {
                return false;
            }
        }
        return true;
    }

    IElementLinear2D * ElementsLinear2D::findElement(const size_t index1, const size_t index2)
    {
        IElementLinear2D * el = nullptr;
        for(const auto & element : m_Elements)
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
        std::vector<double> result(HygroThermFEM::maxNodeIndex(), 0);

#ifdef STL_MULTITHREADING

        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          const auto indexes = aElement->nodeIndexes();
                          const auto vecR = aElement->rightSideVector();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              std::lock_guard guard(mtx);
                              result[indexes[i] - 1] += vecR[i];
                          }
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
        std::vector<NodeFlux> result(HygroThermFEM::maxNodeIndex(), {0, 0});

        std::vector<std::vector<NodeFlux>> fluxes(HygroThermFEM::maxNodeIndex(),
                                                  std::vector<NodeFlux>());

#ifdef STL_MULTITHREADING
        std::mutex mtx;
        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](const auto & aElement) {
                          const auto indexes = aElement->nodeIndexes();
                          const auto flux = aElement->flux();
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              std::lock_guard guard(mtx);
                              fluxes[indexes[i] - 1].push_back(flux[i]);
                          }
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
                                             [&](double lhs, const NodeFlux & a) { return lhs + a.x; })
                             / fluxes[j].size();
            const double y = std::accumulate(fluxes[j].begin(),
                                             fluxes[j].end(),
                                             0.0,
                                             [&](double lhs, const NodeFlux & a) { return lhs + a.y; })
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
