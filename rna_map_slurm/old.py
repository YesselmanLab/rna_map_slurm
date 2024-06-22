def combine_gzipped_fastq(input_files, output_file):
    with gzip.open(output_file, "wb") as output_gz:
        for input_file in input_files:
            with gzip.open(input_file, "rb") as input_gz:
                for line in input_gz:
                    output_gz.write(line)


def process_pair(name_pair, barcode, zip_files, outdir, tmp_dir):
    construct_barcode, pair = name_pair
    count = 0
    for zip_file in zip_files:
        try:
            with zipfile.ZipFile(zip_file, "r") as zip_ref:
                zip_ref.extract(pair[0], f"{tmp_dir}/{count}")
                zip_ref.extract(pair[1], f"{tmp_dir}/{count}")
            count += 1
        except:
            continue
    mate1_files = glob.glob(f"{tmp_dir}/*/{pair[0]}")
    mate2_files = glob.glob(f"{tmp_dir}/*/{pair[1]}")
    if len(mate1_files) == 0:
        log.warning(f"no files found for {construct_barcode}")
        return
    if len(mate1_files) != len(mate2_files):
        log.warning(f"mismatched files for {construct_barcode}")
        return
    log.info(f"{construct_barcode} {len(mate1_files)} {len(mate2_files)}")
    combine_gzipped_fastq(mate1_files, f"{outdir}/{pair[0]}")
    combine_gzipped_fastq(mate2_files, f"{outdir}/{pair[1]}")
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[0]}", shell=True)
    subprocess.call(f"rm -r {tmp_dir}/*/{pair[1]}", shell=True)


def get_file_size(file_path):
    file_path = os.path.realpath(file_path)
    return os.path.getsize(file_path)


@cli.command()
@click.argument("fastq_dir", type=click.Path(exists=True))
@click.option("--output_dir", default=None)
def int_demultiplex_by_split(fastq_dir, output_dir):
    setup_logging(file_name=f"{fastq_dir}/int_demultiplex.log")
    if output_dir is None:
        output_dir = os.getcwd()
    df = pd.read_csv(f"data.csv")
    read_1_path = fastq_dir + "/test_R1.fastq.gz"
    read_2_path = fastq_dir + "/test_R2.fastq.gz"
    barcode_seq = Path(fastq_dir).stem
    df_sub = df.loc[df["barcode_seq"] == barcode_seq]
    if df_sub.empty:
        log.error(f"barcode_seq {barcode_seq} not found in data.csv")
        return
    row = df_sub.iloc[0]
    # get helices from commandline
    helices = []
    args = row["demult_cmd"].split()
    for i in range(0, len(args)):
        if args[i] == "--helix" or args[i] == "-helix":
            helices.append([int(args[i + 1]), int(args[i + 2]), int(args[i + 3])])
    unique_code = random_string(10)
    data_path = f"{output_dir}/{unique_code}"
    df_seq = pd.read_csv(f"inputs/rnas/{row['code']}.csv")
    barcode_demultiplex(
        df_seq, Path(read_2_path), Path(read_1_path), helices, data_path
    )
    zip_path = f"{fastq_dir}/int_demultiplexed.zip"
    flatten_and_zip_directory(data_path, zip_path)
    shutil.rmtree(data_path)


@cli.command()
@click.argument("home_dir", type=click.Path(exists=True))
@click.argument("barcode_seq")
@click.option("--threads", default=1)
@click.option("--tmp_dir", default=None)
def join_int_demultiplex(home_dir, barcode_seq, threads, tmp_dir):
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    # get all the zip files
    dirs = glob.glob(f"{home_dir}/data/split-*")
    zip_files = []
    count = 0
    for d in dirs:
        zip_file = f"{d}/{barcode_seq}/int_demultiplexed.zip"
        if not os.path.isfile(zip_file):
            log.warning(f"{zip_file} does not exist")
            continue
        zip_files.append(zip_file)
    # get the name pairs
    df = pd.read_csv(f"{home_dir}/inputs/data.csv")
    row = df[df["barcode_seq"] == barcode_seq].iloc[0]
    df_barcode = pd.read_json(f"{home_dir}/inputs/barcode_jsons/{row['code']}.json")
    pairs = {}
    for _, r in df_barcode.iterrows():
        full_barcode = r["full_barcode"]
        pairs[full_barcode] = [
            f"{full_barcode}_mate1.fastq.gz",
            f"{full_barcode}_mate2.fastq.gz",
        ]

    outdir = f"{home_dir}/int_demultiplexed/{barcode_seq}"
    os.makedirs(outdir, exist_ok=True)
    tmp_dir = tmp_dir + "/" + random_string(10)
    os.makedirs(tmp_dir, exist_ok=True)

    with ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(
            process_pair,
            pairs.items(),
            [barcode_seq] * len(pairs),
            [zip_files] * len(pairs),
            [outdir] * len(pairs),
            [tmp_dir] * len(pairs),
        )

    count = 0
    for zip_file in zip_files:
        shutil.rmtree(f"{tmp_dir}/{count}")
        count += 1
