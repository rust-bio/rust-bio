use std::env;
use std::path::PathBuf;
use std::path::Path;
use regex::Regex;

pub fn expand_path(path: &Path) -> PathBuf {
    let path_str = path.to_string_lossy();

    // This regex matches occurrences of $VAR or ${VAR}
    let re = Regex::new(r"\$(\w+|\{\w+\})").unwrap();

    let expanded_path = re.replace_all(&path_str, |caps: &regex::Captures| {
        let var = caps[1].trim_matches(|c| c == '{' || c == '}');
        env::var(var).unwrap_or_else(|_| caps[0].to_owned())
    });

    PathBuf::from(expanded_path.into_owned())

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_expand_vars() {
        let path = Path::new("$HOME/test.txt");
        let expanded_path = expand_path(path);
        assert_eq!(expanded_path, PathBuf::from(format!("{}/test.txt", env::var("HOME").unwrap())));
    }
}