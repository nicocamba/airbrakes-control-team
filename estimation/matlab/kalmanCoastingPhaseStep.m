function [input, ctx] = kalmanCoastingPhaseStep(t, x, ctx)
    x = x(:);
    ctx = EstimationCoastingPhase(x, ctx);
    [input, ctx] = kalmanCoastingPhase(t, x, ctx);
    ctx.xprev = x;
end