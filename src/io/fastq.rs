
enum State {
    Name,
    Seq,
    Qual
}

#[deriving(Default)]
struct Record {
    name: String,
    seq: String,
    qual: String
}


fn read_fastq(file) {
    let (tx, rx) = channel();
    spawn(move || {
        let mut record : Record = None;
        for (i, line) in file.lines().enumerate() {
            match i % 4 {
                0 => {
                    if record != None {
                        tx.send(record);
                    }
                    record = Record();
                    record.name = line[1..-1];
                },
                1 => { record.seq = line; },
                2 => { skip; },
                3 => { record.qual = line; }
            }
        }
    });
    rx;
}


fn read_fastq(file) {
    let mut state = State::Name;
    let mut record = Nil;

    for (i, line) in file.lines().enumerate() {
        if line[0] == '>' and state != State::Seq {
            
        }
        match line[0] {
            '>' => {
                state = State::Name;
            },
            '+' => {
                state = State::Qual;
            },
            _ => {
                
            }
        }
        match state {
            State::Name => {
                if line[0] == '>' {
                    record = Record::new(line[1:]);
                    state = State::Seq;
                }
                // TODO error reporting
            },
            State::Seq => {
                if line[0] == '+' {
                    state = State::Qual;
                }
                else {
                    seq += line[:-1];
                }
            },
            State::Qual => {
                
            }
        }
    }
}
