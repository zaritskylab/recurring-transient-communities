from .arcos_wrapper import ArcosParams

PARAMS = ArcosParams(
        # According to ARCOS paper defaults
        smoothK=3,
        binThr=0.4,
        biasMet='runmed',
        peakThr=0.3,
        biasK=51,

        # Small Radius (~2 cell diameters)
        neighborhoodSize=14,
        minClsz=1,
        nPrev=5,

        # No filters
        minDuration=1,
        minTotalEventSize=3
    )