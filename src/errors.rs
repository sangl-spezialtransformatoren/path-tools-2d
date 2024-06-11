#[derive(Debug, Clone)]
pub struct ValueError<'a>(pub &'a str);

#[derive(Debug, Clone)]
pub struct NotImplemented<'a>(pub &'a str);
