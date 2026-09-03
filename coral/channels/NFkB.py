from sieve.channel_set import ChannelSet
from sieve.channels import CompositionBiasChannel
from sieve.channels import ShortMotifChannel
from sieve.channels import NetChargeChannel
from sieve.channels import HydropathyChannel
from sieve.channels import HydrophobicMomentChannel
from sieve.channels import DisorderPropensityChannel
from sieve.channels import SequenceEntropyChannel
from sieve.channels import DipeptideFrequencyChannel

channels = ChannelSet([
    # GRR
    CompositionBiasChannel("GRR_comp", radius=16, residues="G"),
    HydropathyChannel("GRR_hydropathy", radius=16),
    SequenceEntropyChannel("GRR_entropy", radius=16),
    DisorderPropensityChannel("GRR_disorder", radius=16),

    # TAD
    CompositionBiasChannel("TAD_comp", radius=8, residues="FLWV"),
    ShortMotifChannel("TAD_motif", r"[FDAWV]..[ILVWY][FDAWV]..[FDAWV]"),
    DisorderPropensityChannel("TAD_disorder", radius=14),
    NetChargeChannel("TAD_netcharge", radius=14),
    HydrophobicMomentChannel("TAD_hydrophobic_moment", radius=6),  # radius stays 6 for all hydrophobic moment computations

    # NHR
    DipeptideFrequencyChannel("sp_dipeptide", radius=18, dipeptide="SP"),
    ShortMotifChannel("NHR_motif", r"[PV].I.[IT]"),
    HydropathyChannel("NHR_hydropathy", radius=18),
    DisorderPropensityChannel("NHR_disorder", radius=18),
])
