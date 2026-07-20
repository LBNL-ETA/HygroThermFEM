#include <numeric>
#include <stdexcept>

#ifdef STL_MULTITHREADING
#    include <execution>
#endif

#include "Elements2D.hxx"

namespace HygroThermFEM
{
    SquareMatrix ElementsLinear2D::conductanceMatrix(const size_t maxNodeIndex)
    {
        // Assemble into a triplet list and build the sparse matrix once via setFromTriplets (which
        // sums duplicate (row, col) entries -- exactly the FE assembly semantics). Inserting element
        // contributions one at a time with coeffRef into an unreserved sparse matrix is ~O(nnz^2)
        // and dominated the steady-state solve.
        constexpr size_t entriesPerElement{numOfQuadrilateralNodes * numOfQuadrilateralNodes};
        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.reserve(m_Elements.size() * entriesPerElement);

#ifdef STL_MULTITHREADING

        std::mutex mtx;

        std::for_each(std::execution::par_unseq,
                      std::begin(m_Elements),
                      std::end(m_Elements),
                      [&](auto && aElement) {
                          const auto indexes = aElement->nodeIndexes();
                          const auto conductance = aElement->DDuMatrices();
                          const auto condDer = aElement->DpDuMatrices();
                          std::vector<Eigen::Triplet<double>> local;
                          local.reserve(entriesPerElement);
                          for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
                          {
                              for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                              {
                                  local.emplace_back(static_cast<int>(indexes[i] - 1),
                                                     static_cast<int>(indexes[j] - 1),
                                                     conductance(i, j) + condDer(i, j));
                              }
                          }
                          const std::lock_guard<std::mutex> guard{mtx};
                          tripletList.insert(tripletList.end(), local.begin(), local.end());
                      });

#else
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto conductance = element->DDuMatrices();
            const auto condDer = element->DpDuMatrices();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                for(size_t j = 0; j < numOfQuadrilateralNodes; ++j)
                {
                    tripletList.emplace_back(static_cast<int>(indexes[i] - 1),
                                             static_cast<int>(indexes[j] - 1),
                                             conductance(i, j) + condDer(i, j));
                }
            }
        }
#endif

        return SquareMatrix{maxNodeIndex, tripletList};
    }

    std::vector<double> ElementsLinear2D::getLumpedMass(const size_t maxNodeIndex,
                                                        const double DTime) const
    {
        const auto numOfNodes{maxNodeIndex};
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

    SquareMatrix ElementsLinear2D::getMassMatrix(const size_t maxNodeIndex, const double DTime) const
    {
        SquareMatrix Capacitance{maxNodeIndex};

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

    std::vector<double> ElementsLinear2D::RVector(const size_t maxNodeIndex) const
    {
        std::vector<double> result(maxNodeIndex, 0);

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

    void ElementsLinear2D::setVolumetricSource(const double value)
    {
        for(auto & element : m_Elements)
        {
            element->setVolumetricSource(value);
        }
    }

    void ElementsLinear2D::setVolumetricSource(const std::vector<double> & perElement)
    {
        if(perElement.size() != m_Elements.size())
        {
            throw std::runtime_error(
              "Volumetric source vector size does not match the number of elements.");
        }
        for(std::size_t index = 0; index < m_Elements.size(); ++index)
        {
            m_Elements[index]->setVolumetricSource(perElement[index]);
        }
    }

    std::vector<double> ElementsLinear2D::volumetricSourceVector(const size_t maxNodeIndex) const
    {
        std::vector<double> result(maxNodeIndex, 0);
        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto load = element->volumetricSourceVector();
            for(size_t i = 0; i < numOfQuadrilateralNodes; ++i)
            {
                result[indexes[i] - 1] += load[i];
            }
        }
        return result;
    }

    std::vector<NodeFlux> ElementsLinear2D::flux(const size_t maxNodeIndex) const
    {
        std::vector<NodeFlux> result(maxNodeIndex, {0, 0});

        std::vector<std::vector<NodeFlux>> fluxes(maxNodeIndex, std::vector<NodeFlux>());

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
                              if(i < flux.size() && indexes[i] >= 1 && indexes[i] <= maxNodeIndex)
                              {
                                  fluxes[indexes[i] - 1].push_back(flux[i]);
                              }
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
                if(i < flux.size() && indexes[i] >= 1 && indexes[i] <= maxNodeIndex)
                {
                    fluxes[indexes[i] - 1].push_back(flux[i]);
                }
            }
        }
#endif

        // Now need to average them
        for(size_t j = 0; j < fluxes.size(); ++j)
        {
            // A node that belongs to no element (e.g. an auto-enclosure environment node) has no
            // contributing element flux; leave it at zero instead of dividing by zero.
            if(fluxes[j].empty())
            {
                result[j] = {0.0, 0.0};
                continue;
            }
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
