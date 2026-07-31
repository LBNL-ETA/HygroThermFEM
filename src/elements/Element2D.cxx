#include <cassert>
#include <numbers>

#include "Common.hxx"
#include "Element2D.hxx"
#include "IntegrationPoints.hxx"
#include "Nodes.hxx"
#include "Materials.hxx"
#include "QuadrilateralLocal2D.hxx"
#include "VectorOperators.hxx"
#include "SimulationProperties.hxx"

namespace
{
    //! Scales column j of a matrix by factor j, turning the stiffness matrix of a
    //! coefficient into the transport matrix of a product potential (see the two-argument DDu).
    HygroThermFEM::SquareMatrix scaleColumns(const HygroThermFEM::SquareMatrix & matrix,
                                             std::span<const double> factors)
    {
        HygroThermFEM::SquareMatrix scaled{matrix};
        for(std::size_t row = 0; row < scaled.size(); ++row)
        {
            for(std::size_t col = 0; col < scaled.size(); ++col)
            {
                scaled(row, col) *= factors[col];
            }
        }
        return scaled;
    }
}   // namespace

namespace HygroThermFEM
{
    //////////////////////////////////////////////////////////////////////////////
    ///  IQLEMatrix2D
    //////////////////////////////////////////////////////////////////////////////
    IQLEIntegrator2D::IQLEIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        m_Global2D{t_Element},
        m_IntegrationMatrix{numOfQuadrilateralNodes, SquareMatrix{numOfQuadrilateralNodes}}
    {}

    SquareMatrix IQLEIntegrator2D::integrate(std::span<const double> t_Values) const
    {
        const auto count = IntegrationPoints2D::Instance().count2D();

        std::vector<std::vector<double>> aMatrix{numOfQuadrilateralNodes,
                                                 std::vector<double>(numOfQuadrilateralNodes, 0)};
        for(auto integrationPoint = 0u; integrationPoint < count; ++integrationPoint)
        {
            auto & intPointMatrix = m_IntegrationMatrix[integrationPoint];

            for(size_t row = 0; row < t_Values.size(); ++row)
            {
                for(size_t col = 0; col < t_Values.size(); ++col)
                {
                    aMatrix[row][col] += intPointMatrix(row, col)
                                         * 0.5 * (t_Values[row] + t_Values[col]);
                }
            }
        }

        return SquareMatrix{aMatrix};
    }

    SquareMatrix
      IQLEIntegrator2D::integrateInterpolated(std::span<const double> t_Values) const
    {
        const auto count = IntegrationPoints2D::Instance().count2D();
        const auto & localElement = QuadrilateralLinearLocal2D::Instance();

        std::vector<std::vector<double>> aMatrix{numOfQuadrilateralNodes,
                                                 std::vector<double>(numOfQuadrilateralNodes, 0)};
        for(auto integrationPoint = 0u; integrationPoint < count; ++integrationPoint)
        {
            const auto & psi = localElement.Psi(integrationPoint);
            auto coefficient = 0.0;
            for(size_t idx = 0; idx < t_Values.size(); ++idx)
            {
                coefficient += psi[idx] * t_Values[idx];
            }

            auto & intPointMatrix = m_IntegrationMatrix[integrationPoint];
            for(size_t row = 0; row < t_Values.size(); ++row)
            {
                for(size_t col = 0; col < t_Values.size(); ++col)
                {
                    aMatrix[row][col] += intPointMatrix(row, col) * coefficient;
                }
            }
        }

        return SquareMatrix{aMatrix};
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDDUConductance2D
    //////////////////////////////////////////////////////////////////////////////

    QLEDDuIntegrator2D::QLEDDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            const auto & DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto & DPsiDxDyMatrix = m_IntegrationMatrix[integrationPoint];
            for(auto row = 0u; row < DPsiDxDyMatrix.size(); ++row)
            {
                for(auto col = 0u; col < DPsiDxDyMatrix.size(); ++col)
                {
                    DPsiDxDyMatrix(row, col) =
                      (DPsiDx[row] * DPsiDx[col] + DPsiDy[row] * DPsiDy[col]) * det;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLEDpDuIntegrator2D
    //////////////////////////////////////////////////////////////////////////////

    QLEDpDuIntegrator2D::QLEDpDuIntegrator2D(const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {}

    void QLEDpDuIntegrator2D::setIndependentVariables(std::span<const double> t_Values)
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        assert(t_Values.size() == numOfIntegrationPoints);

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & psi = aElement.Psi(integrationPoint);
            const auto & DPsiDx = m_Global2D.DPsiDx(integrationPoint);
            const auto & DPsiDy = m_Global2D.DPsiDy(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            auto gammaX = 0.0;
            auto gammaY = 0.0;
            for(auto idx = 0u; idx < numOfIntegrationPoints; ++idx)
            {
                gammaX += DPsiDx[idx] * t_Values[idx];
                gammaY += DPsiDy[idx] * t_Values[idx];
            }

            std::vector<std::vector<double>> matrix{numOfIntegrationPoints,
                                                    std::vector<double>(numOfIntegrationPoints, 0)};

            for(auto row = 0u; row < numOfIntegrationPoints; ++row)
            {
                for(auto col = 0u; col < numOfIntegrationPoints; ++col)
                {
                    matrix[row][col] =
                      det * (psi[row] * DPsiDx[col] * gammaX + psi[row] * DPsiDy[col] * gammaY);
                }
            }
            m_IntegrationMatrix[integrationPoint] = SquareMatrix{matrix};
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  QLECapacitance2D
    //////////////////////////////////////////////////////////////////////////////

    QLECapacitanceIntegrator2D::QLECapacitanceIntegrator2D(
      const QuadrilateralLinearGlobal2D & t_Element) :
        IQLEIntegrator2D{t_Element}
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        auto & aElement = QuadrilateralLinearLocal2D::Instance();

        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto & psi = aElement.Psi(integrationPoint);
            const auto det = m_Global2D.det(integrationPoint);

            for(auto row = 0u; row < numOfIntegrationPoints; ++row)
            {
                for(auto col = 0u; col < numOfIntegrationPoints; ++col)
                {
                    m_IntegrationMatrix[integrationPoint](row, col) =
                      det * psi[row] * psi[col];
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::IElementLinear2D(Nodes & nodePool,
                                       Materials & materialPool,
                                       const size_t index1,
                                       const size_t index2,
                                       const size_t index3,
                                       const size_t index4,
                                       const std::string & materialName,
                                       const Variable variable,
                                       const bool isLinear) :
        m_Material{materialPool.material(materialName)},
        m_FluxVariable(variable),
        m_Nodes{nodePool.getNode(index1),
                nodePool.getNode(index2),
                nodePool.getNode(index3),
                nodePool.getNode(index4)},
        m_Global2D{nodePool.getNode(index1),
                   nodePool.getNode(index2),
                   nodePool.getNode(index3),
                   nodePool.getNode(index4)},
        m_QLECapacitance2D{m_Global2D},
        m_QLEDDu2D{m_Global2D},
        m_QLEDpDu2D{m_Global2D},
        m_Linear{isLinear && m_Material.isLinear()}
    {
        /// Evaluating material influence in every node (This is important to know when
        /// calculating water content).
        // TODO: This will perform for every domain. Consequence is that weighting factors will be
        // doubled when two subdomains are used or tripled if three domains are used. However,
        // results will be correct since ratio will be correct.
        while(!m_Nodes.last())
        {
            /// Form triangle of nodes. node1 is in center and angle is calculated at that node.
            auto & node1 = nodePool.getNode(m_Nodes.current().getNodeNumber());
            auto & node2 = nodePool.getNode(m_Nodes.previous().getNodeNumber());
            auto & node3 = nodePool.getNode(m_Nodes.next().getNodeNumber());

            /// Weighting coefficient depends on angle that is form by nodes next to node1.
            /// That coefficient is fraction of full circle.
            const auto weightingCoefficient =
              angleBetweenNodes(node1, node2, node3) / (2 * std::numbers::pi);
            /// Node will have possibility to calculate certain properties that will be
            /// material dependent. Pass the already-resolved material reference directly.
            node1.assignMaterial(m_Material, weightingCoefficient);
            m_Nodes.moveToNext();
        }
    }

    SquareMatrix IElementLinear2D::DDuMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};

        for(const auto & term : m_DDuFunctions)
        {
            const auto coefficients = term.coefficient->values(m_Nodes);
            const auto stiffness = m_InterpolateCoefficientsAtGaussPoints
                                     ? m_QLEDDu2D.integrateInterpolated(coefficients)
                                     : m_QLEDDu2D.integrate(coefficients);
            result += term.nodalFactor == nullptr
                        ? stiffness
                        : scaleColumns(stiffness, term.nodalFactor->values(m_Nodes));
        }

        return result;
    }

    SquareMatrix IElementLinear2D::DpDuMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};

        /// The integration matrix has to be rebuilt for every term because the independent
        /// variables change, but setIndependentVariables() replaces it wholesale, so the
        /// integrator itself is reused.
        for(const auto & cond : m_DpDuFunctions)
        {
            const auto aDerivatives = cond.derivativeValue->values(m_Nodes);
            m_QLEDpDu2D.setIndependentVariables(aDerivatives);
            const auto values = cond.fixedValue->values(m_Nodes);
            result += m_InterpolateCoefficientsAtGaussPoints
                        ? m_QLEDpDu2D.integrateInterpolated(values)
                        : m_QLEDpDu2D.integrate(values);
        }

        return result;
    }

    SquareMatrix IElementLinear2D::capacitanceMatrices() const
    {
        SquareMatrix result{numOfQuadrilateralNodes};
        for(const auto & cap : m_CapacitanceFunctions)
        {
            const auto values = cap->values(m_Nodes);
            result += m_LumpCapacityNodally ? nodalLumpedCapacity(values)
                                            : m_QLECapacitance2D.integrate(values);
        }
        for(const auto & cap : m_NodalCapacitanceFunctions)
        {
            result += nodalLumpedCapacity(cap->values(m_Nodes));
        }
        return result;
    }

    SquareMatrix
      IElementLinear2D::nodalLumpedCapacity(std::span<const double> nodalCapacity) const
    {
        // Nodal volumes v_i = integral(psi_i) = row sums of the consistent mass matrix
        // integral(psi_i psi_j) (obtained by integrating with a unit coefficient).
        const auto massMatrix =
          m_QLECapacitance2D.integrate(std::vector<double>(numOfQuadrilateralNodes, 1.0));

        SquareMatrix result{numOfQuadrilateralNodes};
        for(std::size_t i = 0; i < numOfQuadrilateralNodes; ++i)
        {
            double nodalVolume = 0.0;
            for(std::size_t j = 0; j < numOfQuadrilateralNodes; ++j)
            {
                nodalVolume += massMatrix(i, j);
            }
            result(i, i) = nodalVolume * nodalCapacity[i];
        }
        return result;
    }

    std::vector<std::array<double, 2>> IElementLinear2D::stateGradientsAtGaussPoints() const
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        const auto values = m_Nodes.properties(m_FluxVariable);
        assert(values.size() == numOfQuadrilateralNodes);

        std::vector<std::array<double, 2>> gradients(numOfIntegrationPoints, {0.0, 0.0});
        for(size_t gpIdx = 0; gpIdx < numOfIntegrationPoints; ++gpIdx)
        {
            const auto DPsiDx = m_Global2D.DPsiDx(gpIdx);
            const auto DPsiDy = m_Global2D.DPsiDy(gpIdx);
            assert(DPsiDx.size() == numOfIntegrationPoints);
            double DtDx{0};
            double DtDy{0};
            for(size_t nodeIdx = 0; nodeIdx < values.size(); ++nodeIdx)
            {
                DtDx += DPsiDx[nodeIdx] * values[nodeIdx];
                DtDy += DPsiDy[nodeIdx] * values[nodeIdx];
            }
            gradients[gpIdx] = {DtDx, DtDy};
        }
        return gradients;
    }

    std::vector<NodeFlux> IElementLinear2D::flux() const
    {
        const auto gradients = stateGradientsAtGaussPoints();
        std::vector<NodeFlux> results;
        for(const auto & cond : m_ConductanceFunctions)
        {
            const auto conductivityValues = cond->values(m_Nodes);
            assert(conductivityValues.size() == numOfQuadrilateralNodes);

            // Extrapolate flux from 2x2 Gauss points to corner nodes using the
            // standard bilinear extrapolation matrix.  Coefficients derive from
            // evaluating the Gauss-point shape functions at the corner nodes:
            //   A = 1 + sqrt(3)/2,  B = -1/2,  C = 1 - sqrt(3)/2
            constexpr auto gaussToNodeA = 1.866025404;
            constexpr auto gaussToNodeB = -0.5;
            constexpr auto gaussToNodeC = 0.133974596;
            const std::vector<std::vector<double>> extrapolationCoefficients{
              {gaussToNodeA, gaussToNodeB, gaussToNodeC, gaussToNodeB},
              {gaussToNodeB, gaussToNodeA, gaussToNodeB, gaussToNodeC},
              {gaussToNodeC, gaussToNodeB, gaussToNodeA, gaussToNodeB},
              {gaussToNodeB, gaussToNodeC, gaussToNodeB, gaussToNodeA}};

            for(size_t nodeIdx = 0u; nodeIdx < numOfQuadrilateralNodes; ++nodeIdx)
            {
                auto fluxX{0.0};
                auto fluxY{0.0};
                for(const auto coeff : extrapolationCoefficients[nodeIdx])
                {
                    fluxX -= conductivityValues[nodeIdx] * gradients[nodeIdx][0] * coeff;
                    fluxY -= conductivityValues[nodeIdx] * gradients[nodeIdx][1] * coeff;
                }
                results.emplace_back(fluxX, fluxY);
            }
        }
        return results;
    }

    std::vector<NodeFlux> IElementLinear2D::fluxAtGaussPoints() const
    {
        const auto gradients = stateGradientsAtGaussPoints();
        const auto cond = meanConductivity();
        std::vector<NodeFlux> result;
        result.reserve(gradients.size());
        for(const auto & grad : gradients)
        {
            result.emplace_back(-cond * grad[0], -cond * grad[1]);
        }
        return result;
    }

    std::vector<std::array<double, 2>> IElementLinear2D::gaussPointGlobalCoordinates() const
    {
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        std::vector<std::array<double, 2>> coords;
        coords.reserve(numOfIntegrationPoints);
        for(size_t gpIdx = 0; gpIdx < numOfIntegrationPoints; ++gpIdx)
        {
            coords.push_back({m_Global2D.xg(gpIdx), m_Global2D.yg(gpIdx)});
        }
        return coords;
    }

    double IElementLinear2D::meanConductivity() const
    {
        double total{0.0};
        std::size_t count{0};
        for(const auto & cond : m_ConductanceFunctions)
        {
            for(const auto value : cond->values(m_Nodes))
            {
                total += value;
                ++count;
            }
        }
        return count == 0 ? 0.0 : total / static_cast<double>(count);
    }

    INode2D & IElementLinear2D::getNode(const std::size_t index) const
    {
        assert(index < m_Nodes.size());
        return m_Nodes[index];
    }

    std::vector<std::size_t> IElementLinear2D::nodeIndexes() const
    {
        return m_Global2D.nodeIndexes();
    }

    const IMaterial & IElementLinear2D::getMaterial() const
    {
        return m_Material;
    }

    bool IElementLinear2D::haveBothNodes(const size_t index1, const size_t index2) const
    {
        bool node1Found = false;
        bool node2Found = false;
        for(auto & node : m_Nodes)
        {
            node1Found = node1Found || node.get().getNodeNumber() == index1;
            node2Found = node2Found || node.get().getNodeNumber() == index2;
        }
        return node1Found && node2Found;
    }

    double IElementLinear2D::angleBetweenNodes(const Node2D & node1,
                                               const Node2D & node2,
                                               const Node2D & node3) const
    {
        auto angle = std::abs(std::atan2(node3.Y() - node1.Y(), node3.X() - node1.X())
                              - std::atan2(node2.Y() - node1.Y(), node2.X() - node1.X()));
        if(angle > std::numbers::pi)
        {
            angle -= std::numbers::pi;
        }
        return angle;
    }

    bool IElementLinear2D::isLinear() const
    {
        return m_Linear;
    }

    void IElementLinear2D::setVolumetricSource(const double value)
    {
        m_VolumetricSource = value;
    }

    std::vector<double> IElementLinear2D::volumetricSourceVector() const
    {
        // Consistent load vector for a constant source: q * integral(psi_i dA),
        // integrated with the same 2x2 Gauss rule as every other element matrix.
        // Kept separate from rightSideVector so the source enters the steady and
        // transient right hand sides symmetrically without altering what the
        // steady path takes from the elements otherwise.
        std::vector<double> result(numOfQuadrilateralNodes, 0);
        if(m_VolumetricSource == 0.0)
        {
            return result;
        }
        const auto & localElement = QuadrilateralLinearLocal2D::Instance();
        const auto numOfIntegrationPoints = IntegrationPoints2D::Instance().count2D();
        for(std::size_t integrationPoint = 0; integrationPoint < numOfIntegrationPoints;
            ++integrationPoint)
        {
            const auto det = m_Global2D.det(integrationPoint);
            const auto & psi = localElement.Psi(integrationPoint);
            for(std::size_t node = 0; node < numOfQuadrilateralNodes; ++node)
            {
                result[node] += m_VolumetricSource * det * psi[node];
            }
        }
        return result;
    }

    std::vector<double> IElementLinear2D::rightSideVector() const
    {
        std::vector<double> result(numOfQuadrilateralNodes, 0);

        /// SquareMatrix M{numOfQuadrilateralNodes};
        for(const auto & item : m_Matrix_x_Vector)
        {
            /// Calculate functions base on node properties
            const auto values = item.MatrixFunction->values(m_Nodes);
            /// And then integrate them
            auto M = m_InterpolateCoefficientsAtGaussPoints
                       ? m_QLEDDu2D.integrateInterpolated(values)
                       : m_QLEDDu2D.integrate(values);
            auto B = item.VectorFunction ? item.VectorFunction->values(m_Nodes)
                                         : m_Nodes.properties(item.PropertyVector);

            result = result + M * B;
        }

        return result;
    }

    IElementLinear2D::MatrixVector::MatrixVector(iValue && matrixFunction,
                                                 const Variable propertyVector) :
        MatrixFunction(std::move(matrixFunction)), PropertyVector(propertyVector)
    {}

    IElementLinear2D::MatrixVector::MatrixVector(iValue && matrixFunction,
                                                 iValue && vectorFunction) :
        MatrixFunction(std::move(matrixFunction)), VectorFunction(std::move(vectorFunction))
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  IElementLinear2D::DerivativeFunction
    //////////////////////////////////////////////////////////////////////////////

    IElementLinear2D::DerivativeFunction::DerivativeFunction(iValue fixedValue,
                                                             iValue derivativeValue) :
        fixedValue(std::move(fixedValue)), derivativeValue(std::move(derivativeValue))
    {}

    IElementLinear2D::NodalProductFunction::NodalProductFunction(iValue coefficient,
                                                                 iValue nodalFactor) :
        coefficient(std::move(coefficient)), nodalFactor(std::move(nodalFactor))
    {}

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementThermalLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementThermalLinear2D::ElementThermalLinear2D(Nodes & nodePool,
                                                   Materials & materialPool,
                                                   const size_t index1,
                                                   const size_t index2,
                                                   const size_t index3,
                                                   const size_t index4,
                                                   const std::string & materialName) :
        IElementLinear2D(nodePool,
                         materialPool,
                         index1,
                         index2,
                         index3,
                         index4,
                         materialName,
                         Variable::temperature,
                         // The fusion secant capacity depends on the temperature iterate, so
                         // the element must go through the Newton-Raphson path; the linear
                         // shortcut would freeze the secant at its step-start (tangent)
                         // value and drop the latent heat entirely.
                         SimulationProperties::Instance().excludeLatentHeatOfFusion())
    {
        //////////////////////////////////////////////////////////////////////////////////////
        /// Capacitance functions
        //////////////////////////////////////////////////////////////////////////////////////

        if(m_Material.hasDensity())
        {
            const StateValue liquidContent(Variable::liquid);
            const StateValue iceContent(Variable::ice);

            // Vapor content changes are ignored for now
            const auto equivalentDensity = m_Material.density() + liquidContent + iceContent;

            const auto equivalentCapacitance =
              m_Material.heatCapacity() + (iceContent / Constants::Density_Ice) * Constants::Cp_Ice
              + (liquidContent / Constants::Density_Water) * Constants::Cp_Water;

            auto capacitance = equivalentDensity * equivalentCapacitance;

            Cap(capacitance);
        }

        // Latent heat of fusion (D5 Tier 1): the mass-conservative secant of
        // L_f * lambda(T) per kilogram of condensed water, times the condensed content
        // (liquid + ice). Validated against the 1D reference solver's enthalpy-method
        // freezing (Stefan front, freeze--thaw conservation). Replaces the disabled
        // PhaseChange sketch (its correct form previously lived on the deleted
        // IceContentFix branch). The nodes' liquidPercent is rolled from lambda(T) at
        // each accepted timestep (MultiDomain::transient), which feeds the liquid/ice
        // split of the sensible capacitance above. Registered via CapNodal: the fusion
        // secant is a near-delta spike at nodes inside the freezing ramp, and the
        // consistent pairwise-average integration would smear it onto neighbouring
        // rows, booking phantom storage energy there (observed as a stalled Stefan
        // front losing ~85 % of the incoming flux).
        if(m_Material.hasDensity()
           && !SimulationProperties::Instance().excludeLatentHeatOfFusion())
        {
            const StateValue liquidContent(Variable::liquid);
            const StateValue iceContent(Variable::ice);
            CapNodal(FusionSecantCapacity() * (liquidContent + iceContent));
        }

        //////////////////////////////////////////////////////////////////////////
        /// Conductance
        //////////////////////////////////////////////////////////////////////////

        /// material
        // const auto materialConductivity =
        //  TabularFunction1D(m_Material.thermalConductivityMoistureAndTemperatureDependent(),
        //  Variable::water);

        if(SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent()
           && m_Material.hasThermalConductivityMoistureAndTemperatureDependent())
        {
            auto materialConductivity =
              m_Material.thermalConductivityMoistureAndTemperatureDependent();
            DDu(materialConductivity);

            auto matCond{m_Material.thermalConductivityMoistureAndTemperatureDependent()};
            CondFlux(matCond);

            // A conductivity that varies inside the element needs the Gauss-point
            // interpolated coefficient: the pairwise average scales diagonal entries
            // by k_i but off-diagonals by 0.5 (k_i + k_j), so the row sums of the
            // conduction operator no longer vanish and the imbalance acts as a
            // spurious interior heat source wherever k has a gradient (verified on
            // UCRL-ID-106550 problem 2.9a, where it shifts the k(T) solution by ~10%).
            // Same remedy the moisture element applies to its vapour coefficient.
            m_InterpolateCoefficientsAtGaussPoints = true;
        }
        if(!SimulationProperties::Instance().thermalConductivityTemperatureAndMoistureDependent()
           && m_Material.hasThermalConductivityDry())
        {
            auto materialConductivity = Constant(m_Material.thermalConductivityDry());
            DDu(materialConductivity);

            Constant matCond{m_Material.thermalConductivityDry()};
            CondFlux(matCond);
        }


        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (vapor part)
        //////////////////////////////////////////////////////////////////////
        if(m_Material.hasDiffusionResistanceFactor())
        {
            const VaporPermeability delta{m_Material.diffusionResistanceFactor()};
            if(!SimulationProperties::Instance().excludeHeatOfEvaporation())
            {
                // Interior latent term h_lg * div(g_v). Two corrections against the THERMM
                // documents (see ThermSample_StuccoWall for the failure they caused):
                //
                // - SIGN: Theoretical Model eq. 59 defines the enthalpy of vaporization
                //   with a leading MINUS; heatOfEvaporation() returns the positive
                //   magnitude, so the term must be negated here. With the sign dropped,
                //   evaporating regions HEAT instead of cool, which runs away in
                //   near-saturated walls (self-heating far above both boundary
                //   temperatures).
                //
                // - POTENTIAL: the divergence acts on the vapour flux the moisture
                //   equation actually transports, g_v = -delta grad(c_sat * phi). The
                //   porosity-laden stored content (Variable::vapor = c_sat * airPorosity
                //   * phi) is NOT that flux potential: its air porosity grows ~30x when a
                //   region dries, fabricating vapour-content gradients (and thus latent
                //   energy) that the moisture equation never moves.
                //
                // Chain-rule split, grad(c_sat * phi) = phi c_sat'(T) grad T
                //                                       + c_sat grad phi:
                //
                // - The T-part is a genuine conductivity contribution
                //   k_lat = h_lg delta phi c_sat'(T) and must be assembled IMPLICITLY
                //   into the conductance matrix. Kept explicit on the right-hand side it
                //   acts as explicit diffusion in T; in vapour-open materials k_lat
                //   exceeds the dry conductivity (Fiberglass Batts: 0.053 vs 0.043 W/mK
                //   at 20 C), the fixed-point iteration loses contraction on fine meshes
                //   and a mesh-dependent spatial sawtooth grows (ThermSample_StuccoWall
                //   saturated fine-mesh cases).
                //
                // - Only the phi-part remains an explicit right-hand-side coupling to
                //   the moisture field, mirroring the staggered treatment of every
                //   other cross-variable term.
                DDu(HeatOfEvaporation() * delta * StateValue(Variable::humidity)
                    * SaturationDerivative());
                multiplies(Constant(-1) * HeatOfEvaporation() * delta * SaturationFunction(),
                           Variable::humidity);
            }
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conversion from liquid to gas (air part)
        //////////////////////////////////////////////////////////////////////

        /// TODO: Add this later when air pressure equation is added
        // auto waterVaporPressure = SaturationFunction() * StateValue(Variable::humidity);

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from liquid
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeCapillaryConduction()
           && m_Material.hasSorptionCurve() && m_Material.hasLiquidTransportationCurve())
        {
            auto humidity = StateValue(Variable::humidity);
            const TabularDerivativeSmooth sorptionDerivative(m_Material.sorptionCurve(),
                                                             Variable::humidity);
            const LiquidTransportationCurve Dl(m_Material.liquidTransportationCurve(), m_Material);
            auto cd = Constant(-1) * Dl * sorptionDerivative * Constants::Cp_Water;
            DpDu(cd, humidity);
        }

        //////////////////////////////////////////////////////////////////////
        ///  Conduction from vapor
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeVaporDiffusionConduction()
           && m_Material.hasDiffusionResistanceFactor())
        {
            const VaporPermeability delta{m_Material.diffusionResistanceFactor()};
            auto vapCond = Constant(-1) * delta * Constants::Cp_Vapor;
            // Advected enthalpy follows the vapour flux of the moisture equation,
            // g_v = -delta grad(c_sat * phi) -- the flux-consistent potential, not the
            // porosity-laden stored content (see the latent-term note above).
            auto vaporFluxPotential = SaturationFunction() * StateValue(Variable::humidity);

            DpDu(vapCond, vaporFluxPotential);
        }
    }

    //////////////////////////////////////////////////////////////////////////////
    ///  ElementMoistureLinear2D
    //////////////////////////////////////////////////////////////////////////////

    ElementMoistureLinear2D::ElementMoistureLinear2D(Nodes & nodePool,
                                                     Materials & materialPool,
                                                     const size_t index1,
                                                     const size_t index2,
                                                     const size_t index3,
                                                     const size_t index4,
                                                     const std::string & materialName) :
        IElementLinear2D(nodePool, materialPool, index1, index2, index3, index4, materialName, Variable::humidity, false)
    {
        // Lump the moisture capacity nodally (diag(v_i * xi_i)) so the secant capacity
        // conserves stored moisture exactly; the averaged-then-lumped form only conserves
        // approximately and leaks ~1% under strong temperature gradients. See D6-refinement.
        m_LumpCapacityNodally = true;

        // Interpolate transport coefficients at Gauss points instead of the pairwise nodal
        // average: with interpolation the column sums of every transport matrix telescope
        // to exactly zero, so total moisture is conserved for ANY spatially varying
        // coefficient (c_sat(T) under a temperature gradient in particular) -- even at
        // not-fully-converged iterates. The pairwise average was the source of the ~1%
        // closed-strip moisture creation in the D1 gradient diagnostic.
        m_InterpolateCoefficientsAtGaussPoints = true;

        //////////////////////////////////////////////////////////////////////////////
        /// Water vapor diffusion
        //////////////////////////////////////////////////////////////////////////////
        if(m_Material.hasDiffusionResistanceFactor())
        {
            // Non-const: the assembly templates store a copy through unique_ptr<IValue>,
            // which a const-deduced T cannot provide.
            VaporPermeability delta{m_Material.diffusionResistanceFactor()};
            auto conductance = delta * SaturationFunction();

            // The vapour term grad.(delta grad(phi c_sat)) is assembled on its PRODUCT
            // potential: the stiffness matrix of delta with column j scaled by the nodal
            // saturation c_sat,j. Splitting it instead into a moisture-gradient half
            // (delta c_sat grad phi) and a temperature-gradient half (delta phi grad c_sat)
            // integrates two matrices that sum to the same continuous term but do not cancel
            // on the profile where the continuous flux vanishes: a sealed strip then settles
            // near, rather than on, its closed-form equilibrium phi = C / c_sat(T) (second
            // order in the element size, 3.1e-5 on a 20-element strip under 40 -> 20 C).
            // The product form annihilates that profile exactly on any mesh and keeps the
            // exact conservation of the split form. See D1 and
            // tst/units/validation/SealedStrip_SteadyGradient.unit.cxx.
            DDu(delta, SaturationFunction());

            // Function for flux calculations.
            CondFlux(conductance);
        }

        //////////////////////////////////////////////////////////////////////////////
        /// Water liquid transportation
        //////////////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeWaterLiquidTransportation()
           && m_Material.hasSorptionCurve() && m_Material.hasLiquidTransportationCurve())
        {
            auto sorptionDerivative =
              TabularDerivativeSmooth(m_Material.sorptionCurve(), Variable::humidity);
            auto WaterLiquidTransport =
              LiquidTransportationCurve(m_Material.liquidTransportationCurve(), m_Material)
              * sorptionDerivative;
            DDu(WaterLiquidTransport);
        }

        //////////////////////////////////////////////////////////////////////////////
        /// Creating capacitance function
        //////////////////////////////////////////////////////////////////////////////
        if(m_Material.hasSorptionCurve())
        {
            // Mass-conservative (secant) storage capacity: (w(phi) - w(phi_prev)) /
            // (phi - phi_prev). With the lumped capacity matrix this conserves total
            // moisture across a step, unlike the tangent dw/dphi which creates/destroys
            // moisture where the sorption curve is steep (near saturation). See
            // SorptionSecantCapacity and doc/Moisture Governing Equations.md (D6).
            SorptionSecantCapacity secantCapacity(m_Material.sorptionCurve(), Variable::humidity);
            Cap(secantCapacity);
        }

        //////////////////////////////////////////////////////////////////////
        /// Functions for flux calculations
        //////////////////////////////////////////////////////////////////////
        if(!SimulationProperties::Instance().excludeWaterLiquidTransportation()
           && m_Material.hasLiquidTransportationCurve())
        {
            CondFlux(LiquidTransportationCurve(m_Material.liquidTransportationCurve(), m_Material));
        }
    }

}   // namespace HygroThermFEM
