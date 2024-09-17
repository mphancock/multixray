


if __name__ == "__main__":
    pdb_file = Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30/0.pdb")
    cif_file = Path("/wynton/home/sali/mhancock/xray/data/cifs/7mhf/7mhf.cif")
    out_cif_file = Path(Path.home(), "xray/tmp/1.cif")

    f_model, status_array = get_f_model_from_f_obs(
        pdb_file=pdb_file,
        cif_file=cif_file,
        occs=[0.75,0.25],
    )

    # f_obs_array = randomize_amplitude(
    #     f_obs=f_model
    # )

    # # flags_array = f_obs_array.generate_r_free_flags(
    # #     fraction=.1,
    # #     max_free=None
    # # )

    # for i in range(50):
    #     print(f_obs_array.indices()[i], f_obs_array.data()[i], flags_array.data()[i])

    # status_array = get_status_array(
    #     flags_array=flags_array
    # )

    # print(status_array.data())
    # print(flags_array.data())

    write_cif(
        f_obs=f_model,
        status_array=status_array,
        cif_file=out_cif_file
    )

    # for i in range(flags_array.size()):
    #     if flags_array.data()[i]:
    #         status_array.data()[i] = "f"
    #     else:
    #         status_array.data()[i] = "o"


    # # f_model = get_f_model_from_cif(
    # #     pdb_file=pdb_file,
    # #     cif_file=cif_file,
    # #     res=1.0
    # # )

    # # Set flags.
    # # status_array = miller_ops.get_miller_array(
    # #     f_obs_file=cif_file,
    # #     label="_refln.status"
    # # )
    # # flags_array = status_array.customized_copy(data=status_array.data()=="f")
    # # f_obs_array, flags_array = f_model.common_sets(other=flags_array)


