

pub fn trim_newline(s: &mut String) {
    if s.ends_with("\n") {
        s.pop();
    }
}
