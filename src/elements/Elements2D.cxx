#include <algorithm>
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
        //
        // Every element owns a fixed slice of the triplet list, so the loop needs no locking and
        // the triplet order -- and therefore the summation order inside setFromTriplets -- is the
        // same however the work happens to be scheduled.
        constexpr size_t entriesPerElement{numOfQuadrilateralNodes * numOfQuadrilateralNodes};
        std::vector<Eigen::Triplet<double>> tripletList(m_Elements.size() * entriesPerElement);

        const auto assembleElement = [&](const std::unique_ptr<IElementLinear2D> & element) {
            // m_Elements is contiguous, so the element's slot gives its slice offset.
            const auto elementIndex = static_cast<size_t>(&element - m_Elements.data());
            const auto indexes = element->nodeIndexes();
            const auto conductance = element->DDuMatrices();
            const auto condDer = element->DpDuMatrices();

            auto entry = elementIndex * entriesPerElement;
            for(size_t row = 0; row < numOfQuadrilateralNodes; ++row)
            {
                for(size_t col = 0; col < numOfQuadrilateralNodes; ++col)
                {
                    tripletList[entry++] =
                      Eigen::Triplet<double>{static_cast<int>(indexes[row] - 1),
                                             static_cast<int>(indexes[col] - 1),
                                             conductance(row, col) + condDer(row, col)};
                }
            }
        };

#ifdef STL_MULTITHREADING
        // par, not par_unseq: the element kernel allocates, which is not vectorisation-safe.
        std::for_each(
          std::execution::par, std::begin(m_Elements), std::end(m_Elements), assembleElement);
#else
        std::for_each(std::begin(m_Elements), std::end(m_Elements), assembleElement);
#endif

        return SquareMatrix{maxNodeIndex, tripletList};
    }

    std::vector<double> ElementsLinear2D::getLumpedMass(const size_t maxNodeIndex,
                                                        const double DTime) const
    {
        // Row-summed (lumped) capacitance, accumulated straight into the result. This used to
        // scatter into a dense maxNodeIndex x maxNodeIndex scratch matrix and then sum each row,
        // which cost O(maxNodeIndex^2) time and memory to carry the 16 values per element that
        // actually contribute -- at 10k nodes that is ~800 MB and 10^8 additions per call.
        std::vector<double> mass(maxNodeIndex, 0);

        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto capacitance = element->capacitanceMatrices();
            for(size_t row = 0; row < numOfQuadrilateralNodes; ++row)
            {
                for(size_t col = 0; col < numOfQuadrilateralNodes; ++col)
                {
                    mass[indexes[row] - 1] += capacitance(row, col);
                }
            }
        }

        for(auto & nodalMass : mass)
        {
            nodalMass /= DTime;
        }

        return mass;
    }

    SquareMatrix ElementsLinear2D::getMassMatrix(const size_t maxNodeIndex, const double DTime) const
    {
        SquareMatrix Capacitance{maxNodeIndex};

        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto capacitance = element->capacitanceMatrices();
            for(size_t row = 0; row < numOfQuadrilateralNodes; ++row)
            {
                for(size_t col = 0; col < numOfQuadrilateralNodes; ++col)
                {
                    Capacitance(indexes[row] - 1, indexes[col] - 1) +=
                      capacitance(row, col) / DTime;
                }
            }
        }

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

        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto vecR = element->rightSideVector();
            for(size_t idx = 0; idx < numOfQuadrilateralNodes; ++idx)
            {
                result[indexes[idx] - 1] += vecR[idx];
            }
        }

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

        for(const auto & element : m_Elements)
        {
            const auto indexes = element->nodeIndexes();
            const auto elementFlux = element->flux();
            for(size_t idx = 0; idx < numOfQuadrilateralNodes; ++idx)
            {
                if(idx < elementFlux.size() && indexes[idx] >= 1 && indexes[idx] <= maxNodeIndex)
                {
                    fluxes[indexes[idx] - 1].push_back(elementFlux[idx]);
                }
            }
        }

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
