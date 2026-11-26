"""Recheck all listed public theorems and their transitive axiom reports."""
from pathlib import Path
import os
import re
import subprocess
import sys
from tempfile import TemporaryDirectory

root = Path(__file__).resolve().parent
theorems = {
    'MathCommons.lean': {'MathCommons.' + n for n in (
        'inEveryConvexSet', 'allPairings', 'orthogonal', 'determinantValue', 'witnessInHull')},
    'ForestCut.lean': {'MathCommons.ForestCut.' + n for n in (
        'telescoping', 'constantOnInterval', 'singleCutChange', 'totalChangeAtOpportunity',
        'sameEndpointsSameCrossing', 'cutHistoryUnique', 'crossingDeltaInjective', 'crossingDecisionUnique')},
    'ForestTrajectory.lean': {'MathCommons.ForestTrajectory.' + n for n in (
        'swapAtLeft', 'swapAtRight', 'swapCollision', 'swapEmpty', 'swapFull',
        'countSplit', 'restCountSwap', 'swapPreservesNoncrossingCut', 'stepPreservesNoncrossingCut',
        'equalCountSwapCollision', 'nextSubsetUnique', 'trajectoryPreservesCardinality',
        'subsetCountsSingleOpportunity', 'trajectoryUnique', 'edgeSideSeparates',
        'survivingEdgeSameSide', 'forestTrajectoryUnique')},
    'ForestGateSupport.lean': {'MathCommons.ForestGateSupport.' + n for n in (
        'columnChoicesColumns', 'choicePerm', 'emptyMinor', 'rationalTwoByTwo',
        'nonzeroSumTerm', 'nonzeroSignedFactors', 'nonzeroMinorMatching',
        'matchingLength', 'rectangularMinorZero', 'matchingRowSupport', 'matchingColumnSupport', 'twoCoordinateGateSupport',
        'outsideAgreementStep', 'nonzeroMinorStep', 'twoCoordinateNonzeroMinorStep',
        'gaussianGateMinorStep', 'nonzeroMinorCollision', 'nonzeroMinorEmpty', 'nonzeroMinorFull',
        'nonzeroMinorFactorsTrajectory', 'nonzeroMinorHistoriesUnique',
        'nonzeroPathFactors', 'nonzeroAmplitudesUnique')},
    'ForestGateMagnitude.lean': {'MathCommons.ForestGateMagnitude.' + n for n in (
        'normSqZero', 'normSqOne', 'normSqNeg', 'normSqMul', 'normSqSigned', 'unitaryWeights',
        'gaussianOneMul', 'gaussianMulOne', 'gaussianAddZero', 'gaussianZeroAdd', 'singletonMinor',
        'gaussianTwoByTwo', 'scalarSumMapZero', 'minorFirstColumn', 'identityPrefixMinor',
        'lastChoiceMember', 'scalarSumUniqueKey', 'minorLastColumn', 'selectedChoiceMember',
        'minorSelectedColumn', 'blockDetLagrange', 'unitaryDetNorm', 'minorTwoSupportedColumns',
        'decomposeTwoEndpoints', 'filterTwoEndpoints', 'gateSpectatorRow', 'gateEntries',
        'prefixGateMagnitudes', 'gateEndpointSpectatorZero', 'interveningExchangeSign',
        'interveningExchangeNorm', 'signedNot', 'reverseExchangeSign', 'reverseExchangeNorm',
        'gateSpectatorMinor', 'leftStayMinor', 'rightStayMinor', 'oppositeSignsProduct',
        'bothOccupiedMinor', 'orderedGateMagnitudes', 'normSqNonnegative',
        'unitaryProbabilityBounds', 'reverseGate', 'reverseUnitary', 'emptyFilteredMagnitude',
        'fullFilteredMagnitude')},
    'ForestPathProbability.lean': {'MathCommons.ForestPathProbability.' + n for n in (
        'localProbIllegal', 'nonzeroEndpointsLegal', 'stepMagnitudeLeftFirst', 'localProbSwap',
        'stepMagnitudeRightFirst', 'stepMagnitude', 'spectatorsAgreeSpec', 'spectatorsAgreeOfStep',
        'stepMagnitudeGeneral', 'pathMagnitude', 'emptySchedule', 'zeroStepPath', 'collisionStay',
        'nonzeroProbabilityUnique')},
    'ForestPathNormalization.lean': {'MathCommons.ForestPathNormalization.' + n for n in (
        'sumWithAppend', 'sumWithMap', 'sumWithFlatMap', 'sumWithScale', 'sumWithCongr',
        'successorMass', 'stepProbabilityNonnegative', 'successorLegal', 'everyLegalSuccessor',
        'successorsNodup', 'pathProbabilityPrefixCongr', 'extendedPathProbability',
        'historyLegal', 'everyLegalPrefix', 'historyFuture', 'historyDeterminedByPrefix',
        'historiesNodup', 'pathMassOne', 'pathProbabilityNonnegative', 'pathSquaredMagnitudeMassOne')},
    'ForestEndpointMass.lean': {'MathCommons.ForestEndpointMass.' + n for n in (
        'sumWithZero', 'normSumUnique', 'endpointMember', 'endpointHistoriesNodup',
        'endpointNonzeroUnique', 'endpointMagnitude', 'emptyEndpointFiber')},
    'ForestProductMinor.lean': {'MathCommons.ForestProductMinor.' + n for n in (
        'gaussianAddAssoc', 'gaussianAddComm', 'gaussianMulAdd', 'gaussianAddMul',
        'gaussianMulComm', 'gaussianMulAssoc', 'gaussianSumCongr', 'gaussianSumAdd',
        'gaussianSumScale', 'gaussianSumDelta', 'gaussianSumAppend', 'gaussianSumFlatMap',
        'gaussianSumSwap', 'gaussianSumScaleRight', 'matrixMulIdentity', 'matrixIdentityMul',
        'gateProductOne', 'gateMulSpectator', 'gateMulLeft', 'gateMulRight',
        'minorCongrRows', 'minorReplaceAbsent', 'signedAdd', 'signedScale',
        'minorRowAdd', 'minorRowScale', 'minorRowLinearCombination', 'replaceOwnRow',
        'gateMinorNeither', 'gateMinorLeftOnly', 'gateMinorRightOnly', 'productEmptyMinor',
        'productSingletonMinor', 'oneGateProductMagnitude', 'replaceRowsCommute', 'replaceOtherOwnRow',
        'gateAsRowUpdates', 'gateMinorBothExpansion', 'gateMinorBothOfAlternation', 'signedSum',
        'laplaceSum', 'laplaceTailSum', 'minorAsRows', 'rowDetHeadSum',
        'rowDetTailSum', 'productMinorTupleExpansion', 'gaussianSumPairs', 'productMinorTwoRows',
        'cauchyBinetTwo', 'iteratedCauchyBinetTwo', 'gateProductTwoMinorRecurrence', 'cauchyBinetOne')},
    'ForestPairHistory.lean': {'MathCommons.ForestPairHistory.' + n for n in (
        'pairMemberSublist', 'filterSingle', 'pairFilter', 'subsetOfFilter',
        'pairSubsetInjective', 'orderedPairsNodup', 'nodupMapOn', 'pairSubsetsNodup',
        'twoSubsetRepresentative', 'identitySubsetMinor', 'gaussianSumFilter', 'endpointAsSum',
        'pathAmplitudePrefixCongr', 'extendedPathAmplitude', 'successorMinorSum', 'endpointAmplitudeRecurrence',
        'historySizeTwo', 'productPairEndpointAmplitude', 'productPairMagnitude', 'productPairEmptyFiber',
        'sumWithCongrOn', 'sumWithFilter', 'sumWithAdd', 'sumWithSwap',
        'sumWithUniqueKey', 'endpointProbabilityAsSum', 'endpointPairMassOne', 'productPairMassOne')},
}
allowed = {'propext', 'Classical.choice', 'Quot.sound'}
lean = sys.argv[1] if len(sys.argv) == 2 else 'lean'
version = (root / 'lean-toolchain').read_text().strip().split(':v')[1]
output = subprocess.check_output([lean, '--version'], cwd=root, text=True)
if not re.search(r'\bversion ' + re.escape(version) + r'[, )]', output):
    sys.exit(f'Expected Lean {version}, found {output.strip()}')


def check_module(source, expected, build_dir):
    output_file = Path(build_dir, Path(source).stem + '.olean')
    run = subprocess.run([lean, '-DwarningAsError=true', '-o', str(output_file), source],
                         cwd=root, env={**os.environ, 'LEAN_PATH': build_dir},
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(run.stdout, end='')
    if run.returncode:
        sys.exit(run.returncode)
    seen = set()
    for name, deps in re.findall(r"'([\w.]+)' (?:does not depend on any axioms|depends on axioms:\s*\[([^\]]*)\])", run.stdout):
        used = {entry.strip() for entry in deps.split(',') if entry.strip()}
        if used - allowed:
            sys.exit(f'{name} depends on unapproved axioms: {sorted(used - allowed)}')
        seen.add(name)
    if expected - seen:
        sys.exit(f'Missing axiom reports: {sorted(expected - seen)}')


# Later modules import the dependencies built by this same verification run.
with TemporaryDirectory(prefix='math-commons-lean-') as build_dir:
    for source, expected in theorems.items():
        check_module(source, expected, build_dir)
