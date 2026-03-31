open module gammaspike {
    requires beast.base;
    requires beast.fx;
    requires sampled.ancestors;


    exports gammaspike.clockmodel;
    exports gammaspike.distribution;
    exports gammaspike.logger;
    exports gammaspike.operator;
    exports gammaspike.sitemodel;
    exports gammaspike.tree;
    exports gammaspike.util;


    provides beast.base.core.BEASTInterface with
        
        gammaspike.clockmodel.PunctuatedRelaxedClockModel,
        gammaspike.clockmodel.SpikeSize,
        gammaspike.distribution.StumpedTreePrior,
        gammaspike.distribution.BranchSpikePrior,
        gammaspike.distribution.BranchRatePrior,
        gammaspike.distribution.BinomialPrior,
        gammaspike.logger.StumpedTreeLogger,
        gammaspike.logger.TaxonCountLogger,
        gammaspike.logger.SaltativeProportionLogger,
        gammaspike.operator.SpikeFlipOperator,
        gammaspike.operator.SpikeAndSampledAncestorJump,
        gammaspike.operator.SpikeAndSampledAncestorNeighbourJump,
        gammaspike.operator.StubCreator,
        gammaspike.operator.StumpedTreeScaler,
        gammaspike.operator.StubSwapOperator,
        gammaspike.operator.StubBranchOperator,
        gammaspike.operator.StumpedTreeExchange,
        gammaspike.operator.StumpedTreeUniform,
        gammaspike.operator.StumpedTreeConstantDistanceOperator,
        gammaspike.operator.SpikeUpRateDown,
        gammaspike.operator.StubbedLeafToSampledAncestorJump,
        gammaspike.operator.StubTreeStretchOperator,
        gammaspike.operator.StubEpochFlexOperator,
        gammaspike.tree.BirthDeathStubSimulator,
        gammaspike.tree.Stubs,
        gammaspike.sitemodel.ClockMixtureModel,
        gammaspike.logger.SiteCategoryLogger;

}


